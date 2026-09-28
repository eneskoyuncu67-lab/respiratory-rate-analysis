"""
For EVERY subject in bart_edf_correct_order, plots the ENTIRE raw
respiration recording (not just cut-out blocks) with:
  - every trigger marker as a labeled vertical dashed line (so you can
    see exactly what marker text appears where, in context)
  - the 4 resolved blocks shaded in distinct colors, using the NEW
    marker-parsing logic (comma OR period codes, debounced double-taps,
    single-occurrence "ends at next marker" fallback, with a 60s
    plausibility floor on fallback-derived blocks)

This is a VISUAL VERIFICATION tool - it does not change any pipeline
output. One PNG per subject, so you can scroll through all 50 and
confirm the shaded regions actually line up with real blocks, using
the corrected participant numbering/names throughout.
"""
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mne
import numpy as np

from marker_parsing import get_block_segments  # shared with the full pipeline - same overrides included
from config import EDF_DIR, OUT_DIR as BASE_OUT_DIR

OUT_DIR = os.path.join(BASE_OUT_DIR, 'full_recording_marker_plots')
os.makedirs(OUT_DIR, exist_ok=True)

BLOCK_COLORS = {'1_1': '#4C72B0', '1_2': '#DD8452', '2_1': '#55A868', '2_2': '#C44E52'}
BLOCK_LABELS = {'1_1': '1st Certain', '1_2': '1st Uncertain', '2_1': '2nd Certain', '2_2': '2nd Uncertain'}


def get_sorted_markers(raw):
    """Just the raw markers, sorted chronologically, for labeling on the plot."""
    onsets = raw.annotations.onset
    descs = raw.annotations.description
    order = onsets.argsort()
    return onsets[order], [descs[i] for i in order]


files = sorted(glob.glob(os.path.join(EDF_DIR, '*.edf')))
print(f'Plotting {len(files)} recordings...')

for f in files:
    name = os.path.basename(f)
    try:
        raw = mne.io.read_raw_edf(f, preload=True, verbose=False)
        sfreq = raw.info['sfreq']
        sig = raw.get_data(picks=[raw.ch_names.index('Sensor-H:RSP')])[0]
        t = np.arange(len(sig)) / sfreq

        segments = get_block_segments(raw, name)
        onsets, descs = get_sorted_markers(raw)

        fig, ax = plt.subplots(figsize=(18, 5))
        ax.plot(t, sig, 'k-', linewidth=0.4, zorder=1)

        # shade the 4 resolved blocks
        for key, (s, e) in segments.items():
            ax.axvspan(s, e, color=BLOCK_COLORS.get(key, 'gray'), alpha=0.15, zorder=0)
            ax.text((s + e) / 2, ax.get_ylim()[1] if False else np.nanmax(sig),
                    BLOCK_LABELS.get(key, key), ha='center', va='top', fontsize=8,
                    color=BLOCK_COLORS.get(key, 'gray'), fontweight='bold')

        # every marker as a labeled vertical line, in raw chronological order
        for onset, desc in zip(onsets, descs):
            ax.axvline(onset, color='red', linestyle='--', linewidth=0.6, alpha=0.6, zorder=2)
            ax.text(onset, np.nanmin(sig), str(desc), rotation=90, fontsize=6,
                    ha='right', va='bottom', color='red')

        n_resolved = len(segments)
        ax.set_title(f'{name}   ({n_resolved}/4 blocks resolved)', fontsize=11)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Raw RSP')
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f'{os.path.splitext(name)[0]}.png'), dpi=110)
        plt.close(fig)
        print(f'  {name}: {n_resolved}/4 blocks')
    except Exception as e:
        print(f'  {name}: FAILED - {e}')

print(f'\nSaved all plots to: {OUT_DIR}')
