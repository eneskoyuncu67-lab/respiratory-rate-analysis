"""
Interactive manual peak correction - ONLY for blocks that came back NA or
were flagged low-confidence in rsp_rate_by_block.csv (short block or too
few detected breaths). All other blocks are left as computed by
NeuroKit2 (biosppy) - not reviewed here.

For each block in scope: shows the cleaned RSP segment with its current
detected peaks. You correct it by clicking directly on the plot:
  - Left-click NEAR an existing red peak dot  -> removes that peak
  - Left-click ANYWHERE ELSE on the trace      -> adds a peak at the
                                                   nearest local maximum
                                                   within +/-1s of the click
  - Click the "Save & Next" button             -> saves your corrected
                                                   peak list for this block
                                                   and moves to the next one
  - Click "Skip (keep original)"               -> keeps NeuroKit2's
                                                   original result unchanged
Block start/end times are NOT editable here - only individual breath
peaks. Progress is saved incrementally, so closing partway through does
not lose already-corrected blocks.

Run this locally (not via -c/headless) so a real window can open:
    D:\\bart_v1\\neurokit\\python.exe D:\\bart_v1\\manual_correct_peaks.py
"""
import csv
import datetime
import json
import os
import re

import matplotlib
matplotlib.use('TkAgg')  # needs a real display - this script is meant to be run locally, not headless
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import mne
import neurokit2 as nk
from neurokit2.signal import signal_detrend
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

from marker_parsing import get_block_segments  # shared with rsp_rate_by_block.py and manual_correct_single.py

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
IN_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block.csv')
OUT_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block_corrected.csv')

# Stores the ACTUAL corrected peak times (seconds within the block), keyed
# by "filename|label" - shared with manual_correct_single.py. Needed
# because the CSV only stores the final bpm number; without this, revisiting
# an already-corrected block would silently recompute NeuroKit2's original
# (uncorrected) peaks instead of showing what you actually corrected.
PEAKS_STORE_PATH = os.path.join(OUT_DIR, 'manual_peaks_store.json')


def load_peaks_store():
    if os.path.exists(PEAKS_STORE_PATH):
        with open(PEAKS_STORE_PATH, 'r') as f:
            return json.load(f)
    return {}


def save_peaks_store(store):
    with open(PEAKS_STORE_PATH, 'w') as f:
        json.dump(store, f, indent=2)


# Shared, append-only history file with manual_correct_single.py - see that
# script's comment for why this exists (answers "what did I last work on").
LOG_PATH = os.path.join(OUT_DIR, 'correction_log.csv')
LOG_FIELDS = ['timestamp', 'tool', 'file', 'label', 'action', 'rate_bpm',
              'n_original', 'n_added', 'n_removed', 'edit_ratio']


def append_log(row):
    is_new = not os.path.exists(LOG_PATH)
    with open(LOG_PATH, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=LOG_FIELDS)
        if is_new:
            w.writeheader()
        row = dict(row)
        row['timestamp'] = datetime.datetime.now().isoformat(timespec='seconds')
        row.setdefault('tool', 'batch')
        w.writerow(row)


def print_last_log_entries(n=5):
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
    greedy nearest-neighbor matching within MATCH_TOLERANCE_S. Returns
    edit_ratio = (n_added + n_removed) / n_original - the fraction of the
    algorithm's original output you had to touch by hand.
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
    edit_ratio = (n_added + n_removed) / n_original if n_original > 0 else (1.0 if n_added > 0 else 0.0)

    return {'n_original': n_original, 'n_final': n_final,
            'n_added': n_added, 'n_removed': n_removed, 'edit_ratio': edit_ratio}


CODE_RE = re.compile(r'^([12]),([12])$')
LABEL_TO_KEY = {'1st Uncertain': '1_2', '1st Certain': '1_1',
                '2nd Uncertain': '2_2', '2nd Certain': '2_1'}
LABELS = list(LABEL_TO_KEY.keys())

CLICK_TOLERANCE_S = 1.0   # click within this many seconds of a peak removes it
ADD_SEARCH_WINDOW_S = 1.0  # when adding, search for the local max within +/- this


# get_block_segments is now imported from marker_parsing.py (see top of
# file) - shared with rsp_rate_by_block.py and manual_correct_single.py so
# all three scripts always agree on block boundaries, including the
# hand-diagnosed per-subject overrides.


# ---------------------------------------------------------------------------
# Build the review queue: only NA or low-confidence blocks
# ---------------------------------------------------------------------------
df = pd.read_csv(IN_CSV)

# resume support: if a corrected CSV already exists, start from it instead
if os.path.exists(OUT_CSV):
    df = pd.read_csv(OUT_CSV)
    print(f'Resuming from existing {OUT_CSV}')

if 'corrected' not in df.columns:
    for label in LABELS:
        df[f'{label}_corrected'] = False

HIGH_RATE_THRESHOLD = 18  # also review any block at or above this, even if it
                          # wasn't NA/low-confidence - per explicit request

queue = []
for idx, r in df.iterrows():
    for label in LABELS:
        val = r[label]
        needs_review = (
            pd.isna(val)
            or bool(r.get(f'{label}_low_conf', False))
            or (pd.notna(val) and val > HIGH_RATE_THRESHOLD)
        )
        already_done = bool(r.get(f'{label}_corrected', False))
        if needs_review and not already_done:
            queue.append((idx, r['file'], label))

print(f'{len(queue)} blocks to review this pass '
      f'(NA + low-confidence + >{HIGH_RATE_THRESHOLD} bpm, excluding already-corrected).')
print_last_log_entries()

if not queue:
    print('Nothing left to review.')
    raise SystemExit


# ---------------------------------------------------------------------------
# Interactive session state
# ---------------------------------------------------------------------------
state = {'q_idx': 0, 'peaks_s': None, 't': None, 'sig': None, 'sfreq': None,
         'row_idx': None, 'label': None, 'original_peaks_s': None}

fig, ax = plt.subplots(figsize=(12, 5))
plt.subplots_adjust(bottom=0.2)
ax_save = plt.axes([0.72, 0.05, 0.12, 0.075])
ax_skip = plt.axes([0.85, 0.05, 0.12, 0.075])
btn_save = Button(ax_save, 'Save & Next')
btn_skip = Button(ax_skip, 'Skip (keep original)')

_raw_cache = {}


def load_block(row_idx, name, label):
    if name not in _raw_cache:
        raw = mne.io.read_raw_edf(os.path.join(EDF_DIR, name), preload=True, verbose=False)
        sfreq = raw.info['sfreq']
        sig = raw.get_data(picks=[raw.ch_names.index('Sensor-H:RSP')])[0]
        segments = get_block_segments(raw, name)
        _raw_cache[name] = (sig, sfreq, segments)
    sig, sfreq, segments = _raw_cache[name]
    key = LABEL_TO_KEY[label]
    if key not in segments:
        return None
    start_s, end_s = segments[key]
    i0, i1 = int(start_s * sfreq), int(end_s * sfreq)
    seg = sig[i0:i1]
    if len(seg) < 10:
        return sfreq, seg, [], []
    cleaned = clean_rsp(seg, sfreq)

    # ALWAYS compute NeuroKit2's original automatic detection - this is the
    # fixed baseline used to measure edit_ratio, regardless of whether a
    # manual correction already exists for this block.
    t = np.arange(len(cleaned)) / sfreq
    try:
        _, info = nk.rsp_peaks(cleaned, sampling_rate=sfreq, method='biosppy')
        original_peaks_s = list(t[info['RSP_Peaks']])
    except Exception:
        original_peaks_s = []

    store = load_peaks_store()
    store_key = f'{name}|{label}'
    if store_key in store:
        # already manually corrected before (in this tool or the single-block
        # tool) - load that instead of showing the original detection
        peaks_s = list(store[store_key])
        return sfreq, cleaned, peaks_s, original_peaks_s

    peaks_s = list(original_peaks_s)
    return sfreq, cleaned, peaks_s, original_peaks_s


def draw(reload_block=True):
    """
    reload_block=True: moving to a new queue item - (re)load the signal
        and reset peaks_s to NeuroKit2's fresh detection for that block.
    reload_block=False: just re-rendering the SAME block after the user
        added/removed a peak - must NOT touch state['peaks_s'], or every
        manual edit gets silently wiped out on the next redraw (this was
        the bug: draw() used to always reload, so clicks never stuck).
    """
    ax.clear()
    q_idx = state['q_idx']
    if q_idx >= len(queue):
        ax.text(0.5, 0.5, 'All blocks reviewed. Close this window.',
                ha='center', va='center', fontsize=14, transform=ax.transAxes)
        fig.canvas.draw_idle()
        return

    row_idx, name, label = queue[q_idx]

    if reload_block:
        result = load_block(row_idx, name, label)
        if result is None:
            # no marker pair at all for this block - nothing to correct, auto-skip
            advance_and_mark(row_idx, label, np.nan, corrected=False)
            return
        sfreq, cleaned, peaks_s, original_peaks_s = result
        t = np.arange(len(cleaned)) / sfreq
        state.update(row_idx=row_idx, label=label, t=t, sig=cleaned, sfreq=sfreq,
                     peaks_s=list(peaks_s), original_peaks_s=original_peaks_s)

    t = state['t']
    cleaned = state['sig']
    ax.plot(t, cleaned, 'k-', linewidth=0.7)
    if state['peaks_s']:
        ys = np.interp(state['peaks_s'], t, cleaned)
        ax.scatter(state['peaks_s'], ys, color='red', s=40, zorder=3)
    n_breaths = len(state['peaks_s'])
    dur_min = (t[-1] if len(t) else 0) / 60
    rate = 60 * n_breaths / (t[-1]) if len(t) and t[-1] > 0 and n_breaths >= 2 else float('nan')
    ax.set_title(f'[{q_idx+1}/{len(queue)}] {name}  |  {label}  |  '
                 f'{n_breaths} breaths in {dur_min:.2f} min  ->  {rate:.1f} bpm\n'
                 f'Click near a red dot to remove it, click elsewhere to add a peak',
                 fontsize=10)
    ax.set_xlabel('Time in block (s)')
    ax.set_ylabel('Cleaned RSP')
    fig.canvas.draw_idle()


def advance_and_mark(row_idx, label, rate, corrected, stats=None):
    df.at[row_idx, label] = rate
    df.at[row_idx, f'{label}_corrected'] = corrected
    df.at[row_idx, f'{label}_low_conf'] = False if corrected else df.at[row_idx, f'{label}_low_conf']
    if stats is not None:
        df.at[row_idx, f'{label}_n_original'] = stats['n_original']
        df.at[row_idx, f'{label}_n_added'] = stats['n_added']
        df.at[row_idx, f'{label}_n_removed'] = stats['n_removed']
        df.at[row_idx, f'{label}_edit_ratio'] = round(stats['edit_ratio'], 4)
    df.to_csv(OUT_CSV, index=False)
    state['q_idx'] += 1
    draw()


def on_click(event):
    if event.inaxes != ax or state['t'] is None:
        return
    click_t = event.xdata
    if click_t is None:
        return
    peaks_s = state['peaks_s']

    # is the click near an existing peak? -> remove it
    if peaks_s:
        dists = [abs(p - click_t) for p in peaks_s]
        min_i = int(np.argmin(dists))
        if dists[min_i] <= CLICK_TOLERANCE_S:
            peaks_s.pop(min_i)
            draw(reload_block=False)
            return

    # otherwise -> add a peak at the local max within the search window
    t = state['t']
    sig = state['sig']
    lo = max(0, click_t - ADD_SEARCH_WINDOW_S)
    hi = click_t + ADD_SEARCH_WINDOW_S
    mask = (t >= lo) & (t <= hi)
    if not np.any(mask):
        return
    local_t = t[mask]
    local_sig = sig[mask]
    new_peak_t = local_t[np.argmax(local_sig)]
    peaks_s.append(float(new_peak_t))
    peaks_s.sort()
    draw(reload_block=False)


def on_save(event):
    row_idx = state['row_idx']
    label = state['label']
    peaks_s = state['peaks_s']
    if row_idx is None:
        return
    if len(peaks_s) >= 2:
        rate = 60 / np.mean(np.diff(peaks_s))
    else:
        rate = float('nan')

    stats = compute_edit_stats(state['original_peaks_s'], peaks_s)

    # persist the actual peak positions, not just the derived rate - this
    # is what makes revisiting this block later (in either tool) show your
    # correction instead of NeuroKit2's original detection
    _, name, _ = queue[state['q_idx']]
    store = load_peaks_store()
    store[f'{name}|{label}'] = peaks_s
    save_peaks_store(store)

    advance_and_mark(row_idx, label, rate, corrected=True, stats=stats)
    append_log({'tool': 'batch', 'file': name, 'label': label, 'action': 'saved',
                'rate_bpm': round(rate, 2) if not np.isnan(rate) else '',
                'n_original': stats['n_original'], 'n_added': stats['n_added'],
                'n_removed': stats['n_removed'], 'edit_ratio': round(stats['edit_ratio'], 4)})


def on_skip(event):
    row_idx = state['row_idx']
    label = state['label']
    if row_idx is None:
        state['q_idx'] += 1
        draw()
        return
    original_rate = df.at[row_idx, label]
    # "skip" means keep the original algorithmic result untouched -> 0 edits
    n_orig = len(state['original_peaks_s']) if state['original_peaks_s'] is not None else 0
    stats = {'n_original': n_orig, 'n_final': n_orig, 'n_added': 0, 'n_removed': 0, 'edit_ratio': 0.0}
    _, skip_name, _ = queue[state['q_idx']]
    advance_and_mark(row_idx, label, original_rate, corrected=False, stats=stats)
    append_log({'tool': 'batch', 'file': skip_name, 'label': label, 'action': 'skipped'})


fig.canvas.mpl_connect('button_press_event', on_click)
btn_save.on_clicked(on_save)
btn_skip.on_clicked(on_skip)

draw()
plt.show()

print(f'\nSaved progress to: {OUT_CSV}')
