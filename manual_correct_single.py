"""
Interactive peak correction for ONE subject/block at a time, chosen by
you typing it in - not an automatic queue. Loops: enter a subject, enter
a block code, correct peaks in a popup window, save, repeat.

Block code convention (block-first, matches the marker convention used
throughout this project): "1.1"=1st Certain, "1.2"=1st Uncertain,
"2.1"=2nd Certain, "2.2"=2nd Uncertain.

Run this locally so a real window can open:
    D:\\bart_v1\\neurokit\\python.exe D:\\bart_v1\\manual_correct_single.py

At each prompt:
  Subject: type the leading row number (e.g. "13") or part of the name
            (e.g. "erenali") - whichever matches one file uniquely.
  Block:   1.1 / 1.2 / 2.1 / 2.2
  (leave subject blank and press Enter to quit)

In the popup window:
  - Left-click NEAR a red dot  -> removes that peak
  - Left-click ANYWHERE ELSE   -> adds a peak at the nearest local max
                                   within +/-1s of the click
  - "Save & Close"             -> writes the corrected rate back into the
                                   CSV and closes the window so you can
                                   pick the next subject/block
  - "Cancel (discard)"         -> closes without saving anything
"""
import csv
import datetime
import glob
import json
import os
import re

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import mne
import neurokit2 as nk
from neurokit2.signal import signal_detrend
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

from marker_parsing import get_block_segments  # shared with rsp_rate_by_block.py and manual_correct_peaks.py

# Bandpass filter cutoffs - same evidence-based cutoffs as rsp_rate_by_block.py:
# Ruckert-Eheberg et al. (2025, PLOS ONE), 5th/95th percentile respiratory
# rate (12.06 / 20.06 breaths/min) from 2,224 general-population adults
# (KORA-FF4 study), converted to Hz (bpm/60).
RSP_LOWCUT_HZ = 12.06 / 60
RSP_HIGHCUT_HZ = 20.06 / 60
FILTER_ORDER = 2  # matches BioSPPy's own filter order (Carreiras et al., 2015)


def clean_rsp(seg, sfreq):
    """Same filter architecture as NeuroKit2's biosppy cleaning method
    (2nd-order Butterworth bandpass + constant detrend), but with our own
    evidence-based cutoffs instead of BioSPPy's hardcoded 0.1-0.35 Hz."""
    freq = 2 * np.array([RSP_LOWCUT_HZ, RSP_HIGHCUT_HZ]) / sfreq
    b, a = butter(N=FILTER_ORDER, Wn=freq, btype='bandpass')
    filtered = filtfilt(b, a, seg)
    return signal_detrend(filtered, order=0)


from config import EDF_DIR, OUT_DIR
BASE_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block.csv')
OUT_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block_corrected.csv')

# Stores the ACTUAL corrected peak times (seconds within the block), keyed
# by "filename|label". This is what was missing before: the CSV only ever
# stored the final bpm number, not the peak positions themselves, so
# reopening an already-corrected block had nothing to load and silently
# recomputed NeuroKit2's original (uncorrected) detection from scratch.
# Shared with manual_correct_peaks.py so corrections made in either tool
# are visible in both.
PEAKS_STORE_PATH = os.path.join(OUT_DIR, 'manual_peaks_store.json')


def load_peaks_store():
    if os.path.exists(PEAKS_STORE_PATH):
        with open(PEAKS_STORE_PATH, 'r') as f:
            return json.load(f)
    return {}


def save_peaks_store(store):
    with open(PEAKS_STORE_PATH, 'w') as f:
        json.dump(store, f, indent=2)


# A running, append-only history of every save/discard action, in
# chronological order - this answers "which participant/block did I last
# work on", since neither the CSV nor the peaks store keep any history
# (they only ever hold the current end-state). Never overwritten, only
# appended to, so it survives across every session of either tool.
LOG_PATH = os.path.join(OUT_DIR, 'correction_log.csv')
LOG_FIELDS = ['timestamp', 'tool', 'file', 'label', 'action', 'rate_bpm',
              'n_original', 'n_added', 'n_removed', 'edit_ratio']


def append_log(row):
    """row is a dict with (a subset of) LOG_FIELDS as keys - missing keys are left blank."""
    is_new = not os.path.exists(LOG_PATH)
    with open(LOG_PATH, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=LOG_FIELDS)
        if is_new:
            w.writeheader()
        row = dict(row)
        row['timestamp'] = datetime.datetime.now().isoformat(timespec='seconds')
        row.setdefault('tool', 'single')
        w.writerow(row)


def print_last_log_entries(n=5):
    """Shown at startup so you immediately see where you left off last time."""
    if not os.path.exists(LOG_PATH):
        print('No correction history yet - this will be your first entry.')
        return
    with open(LOG_PATH, 'r', newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print('No correction history yet - this will be your first entry.')
        return
    print(f'Last {min(n, len(rows))} action(s) from your correction history:')
    for r in rows[-n:]:
        print(f'  {r["timestamp"]}  [{r["tool"]}]  {r["file"]:28s} {r["label"]:15s} '
              f'{r["action"]:9s} rate={r["rate_bpm"]}')
    print()


MATCH_TOLERANCE_S = 0.5  # a final peak within this of an original peak counts as "kept", not added/removed


def compute_edit_stats(original_peaks_s, final_peaks_s):
    """
    Compares the final (possibly manually-edited) peak list against
    NeuroKit2's original automatic detection for the same block, via
    greedy nearest-neighbor matching within MATCH_TOLERANCE_S.

    Returns counts of how many peaks were removed (present in the
    original detection, absent from the final list) and added (present
    in the final list, not near anything in the original detection),
    plus edit_ratio = (n_added + n_removed) / n_original - i.e. what
    fraction of the algorithm's original output you had to touch by
    hand. If the algorithm originally found zero peaks, edit_ratio is 1.0
    if you added any (100% manual) or 0.0 if you left it empty.
    """
    orig = sorted(original_peaks_s)
    final = sorted(final_peaks_s)
    matched_orig = set()
    matched_final = set()
    for i, o in enumerate(orig):
        best_j, best_d = None, MATCH_TOLERANCE_S
        for j, f_ in enumerate(final):
            if j in matched_final:
                continue
            d = abs(o - f_)
            if d <= best_d:
                best_d, best_j = d, j
        if best_j is not None:
            matched_orig.add(i)
            matched_final.add(best_j)

    n_original = len(orig)
    n_final = len(final)
    n_removed = n_original - len(matched_orig)
    n_added = n_final - len(matched_final)

    if n_original > 0:
        edit_ratio = (n_added + n_removed) / n_original
    else:
        edit_ratio = 1.0 if n_added > 0 else 0.0

    return {'n_original': n_original, 'n_final': n_final,
            'n_added': n_added, 'n_removed': n_removed, 'edit_ratio': edit_ratio}

CODE_RE = re.compile(r'^([12]),([12])$')
CODE_TO_KEY = {'1.1': '1_1', '1.2': '1_2', '2.1': '2_1', '2.2': '2_2'}
KEY_TO_LABEL = {'1_2': '1st Uncertain', '1_1': '1st Certain',
                '2_2': '2nd Uncertain', '2_1': '2nd Certain'}

CLICK_TOLERANCE_S = 1.0
ADD_SEARCH_WINDOW_S = 1.0


# get_block_segments is now imported from marker_parsing.py (see top of
# file) - shared with rsp_rate_by_block.py and manual_correct_peaks.py so
# all three scripts always agree on block boundaries, including the
# hand-diagnosed per-subject overrides.


def find_file(query):
    files = sorted(glob.glob(os.path.join(EDF_DIR, '*.edf')))
    query = query.strip().lower()
    matches = [f for f in files if os.path.basename(f).lower().startswith(query)
               or query in os.path.basename(f).lower()]
    if len(matches) == 0:
        print(f'  No file matches "{query}".')
        return None
    if len(matches) > 1:
        print(f'  "{query}" matches multiple files, be more specific:')
        for m in matches:
            print('   -', os.path.basename(m))
        return None
    return matches[0]


def load_df():
    if os.path.exists(OUT_CSV):
        return pd.read_csv(OUT_CSV)
    df = pd.read_csv(BASE_CSV)
    for label in KEY_TO_LABEL.values():
        df[f'{label}_corrected'] = False
    return df


def correct_one(edf_path, label):
    name = os.path.basename(edf_path)
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)
    sfreq = raw.info['sfreq']
    sig = raw.get_data(picks=[raw.ch_names.index('Sensor-H:RSP')])[0]
    segments = get_block_segments(raw, name)
    key = [k for k, v in KEY_TO_LABEL.items() if v == label][0]

    if key not in segments:
        print(f'  {name} has no marker pair for {label} - nothing to show.')
        return

    start_s, end_s = segments[key]
    i0, i1 = int(start_s * sfreq), int(end_s * sfreq)
    seg = sig[i0:i1]
    if len(seg) < 10:
        print(f'  {name} {label} segment is degenerate ({len(seg)} samples) - nothing to show.')
        return

    cleaned = clean_rsp(seg, sfreq)
    t = np.arange(len(cleaned)) / sfreq

    # ALWAYS compute NeuroKit2's original automatic detection, even if a
    # manual correction already exists - this is the fixed baseline used
    # to measure how much you had to edit (edit_ratio), regardless of how
    # many editing sessions this block has already been through.
    try:
        _, info = nk.rsp_peaks(cleaned, sampling_rate=sfreq, method='biosppy')
        original_peaks_s = list(t[info['RSP_Peaks']])
    except Exception:
        original_peaks_s = []

    store = load_peaks_store()
    store_key = f'{name}|{label}'
    if store_key in store:
        # A previous correction exists for this exact block - load it
        # instead of showing the original detection, so your earlier
        # edits are actually shown, not silently discarded.
        peaks_s = list(store[store_key])
        print(f'  Loaded previously corrected peaks ({len(peaks_s)} breaths).')
    else:
        peaks_s = list(original_peaks_s)

    fig, ax = plt.subplots(figsize=(12, 5))
    plt.subplots_adjust(bottom=0.2)
    ax_save = plt.axes([0.72, 0.05, 0.12, 0.075])
    ax_cancel = plt.axes([0.85, 0.05, 0.12, 0.075])
    btn_save = Button(ax_save, 'Save & Close')
    btn_cancel = Button(ax_cancel, 'Cancel (discard)')
    result = {'saved': False}

    def render():
        ax.clear()
        ax.plot(t, cleaned, 'k-', linewidth=0.7)
        if peaks_s:
            ys = np.interp(peaks_s, t, cleaned)
            ax.scatter(peaks_s, ys, color='red', s=40, zorder=3)
        n_breaths = len(peaks_s)
        dur_min = t[-1] / 60 if len(t) else 0
        rate = 60 / np.mean(np.diff(peaks_s)) if n_breaths >= 2 else float('nan')
        ax.set_title(f'{name}  |  {label}  |  {n_breaths} breaths in {dur_min:.2f} min  '
                     f'->  {rate:.1f} bpm\n'
                     f'Click near a red dot to remove it, click elsewhere to add a peak',
                     fontsize=10)
        ax.set_xlabel('Time in block (s)')
        ax.set_ylabel('Cleaned RSP')
        fig.canvas.draw_idle()

    def on_click(event):
        if event.inaxes != ax or event.xdata is None:
            return
        click_t = event.xdata
        if peaks_s:
            dists = [abs(p - click_t) for p in peaks_s]
            min_i = int(np.argmin(dists))
            if dists[min_i] <= CLICK_TOLERANCE_S:
                peaks_s.pop(min_i)
                render()
                return
        lo, hi = max(0, click_t - ADD_SEARCH_WINDOW_S), click_t + ADD_SEARCH_WINDOW_S
        mask = (t >= lo) & (t <= hi)
        if not np.any(mask):
            return
        new_peak_t = t[mask][np.argmax(cleaned[mask])]
        peaks_s.append(float(new_peak_t))
        peaks_s.sort()
        render()

    def on_save(event):
        result['saved'] = True
        plt.close(fig)

    def on_cancel(event):
        plt.close(fig)

    fig.canvas.mpl_connect('button_press_event', on_click)
    btn_save.on_clicked(on_save)
    btn_cancel.on_clicked(on_cancel)
    render()
    plt.show()

    if result['saved']:
        rate = 60 / np.mean(np.diff(peaks_s)) if len(peaks_s) >= 2 else float('nan')
        stats = compute_edit_stats(original_peaks_s, peaks_s)

        # persist the actual peak positions - this is what makes re-opening
        # this block later show your correction instead of the original
        store = load_peaks_store()
        store[store_key] = peaks_s
        save_peaks_store(store)

        df = load_df()
        row_mask = df['file'] == name
        if not row_mask.any():
            print(f'  WARNING: {name} not found in {OUT_CSV} - rate not saved to CSV '
                  f'(but peaks were saved to {PEAKS_STORE_PATH}).')
            return
        idx = df.index[row_mask][0]
        df.at[idx, label] = rate
        df.at[idx, f'{label}_corrected'] = True
        if f'{label}_low_conf' in df.columns:
            df.at[idx, f'{label}_low_conf'] = False
        df.at[idx, f'{label}_n_original'] = stats['n_original']
        df.at[idx, f'{label}_n_added'] = stats['n_added']
        df.at[idx, f'{label}_n_removed'] = stats['n_removed']
        df.at[idx, f'{label}_edit_ratio'] = round(stats['edit_ratio'], 4)
        df.to_csv(OUT_CSV, index=False)
        print(f'  Saved: {name} {label} -> {rate:.1f} bpm ({len(peaks_s)} breaths)  '
              f'[edited {stats["n_added"]} added / {stats["n_removed"]} removed '
              f'of {stats["n_original"]} original -> ratio {stats["edit_ratio"]:.2f}]')
        append_log({'tool': 'single', 'file': name, 'label': label, 'action': 'saved',
                    'rate_bpm': round(rate, 2) if not np.isnan(rate) else '',
                    'n_original': stats['n_original'], 'n_added': stats['n_added'],
                    'n_removed': stats['n_removed'], 'edit_ratio': round(stats['edit_ratio'], 4)})
    else:
        print('  Discarded - no change saved.')
        append_log({'tool': 'single', 'file': name, 'label': label, 'action': 'discarded'})


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
print('Manual single-block peak correction. Leave subject blank to quit.\n')
print_last_log_entries()
while True:
    subj = input('Subject (row number or name substring): ').strip()
    if not subj:
        break
    edf_path = find_file(subj)
    if edf_path is None:
        continue

    block = input('Block (1.1 / 1.2 / 2.1 / 2.2): ').strip()
    if block not in CODE_TO_KEY:
        print('  Not a valid block code, try again.')
        continue
    label = KEY_TO_LABEL[CODE_TO_KEY[block]]

    correct_one(edf_path, label)
    print()

print('Done.')
