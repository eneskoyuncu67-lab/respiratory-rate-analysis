"""
Plot a raw signal channel from a BART EDF recording with trigger markers overlaid.
No cleaning, no processing - just the raw trace and the markers, for visual inspection.

Usage:
    python plot_with_markers.py <path_to_edf> [channel_name]

Example:
    python plot_with_markers.py "data/edf/41_subjectname_2.edf" "Sensor-B:EEG"

If channel_name is omitted, defaults to Sensor-B:EEG (the cardiac channel - device
mislabels it EEG, but it is actually ECG based on visual QRS morphology).
"""

import sys

import matplotlib.pyplot as plt
import mne
import numpy as np


def plot_with_markers(edf_path, channel='Sensor-B:EEG'):
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)
    sfreq = raw.info['sfreq']
    idx = raw.ch_names.index(channel)
    sig = raw.get_data(picks=[idx])[0]
    t = np.arange(len(sig)) / sfreq

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.plot(t, sig, linewidth=0.4)
    for onset, desc in zip(raw.annotations.onset, raw.annotations.description):
        ax.axvline(onset, color='red', linestyle='--', alpha=0.7)
        ax.text(onset, ax.get_ylim()[1], desc, rotation=90, va='top', fontsize=8, color='red')
    ax.set_xlabel('time (s)')
    ax.set_title(f'{edf_path} - {channel} - raw, unprocessed, with markers')
    plt.tight_layout()
    return fig


if __name__ == '__main__':
    edf_path = sys.argv[1]
    channel = sys.argv[2] if len(sys.argv) > 2 else 'Sensor-B:EEG'
    fig = plot_with_markers(edf_path, channel)
    plt.show()
