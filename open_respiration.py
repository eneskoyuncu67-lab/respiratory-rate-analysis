"""
Open one participant's respiration (RSP) signal and visualize it with NeuroKit2.
No cleaning, no peak detection, no rate computation - just load and look.

Usage:
    python open_respiration.py <path_to_edf>

Example:
    python open_respiration.py "data/edf/01_subjectname.edf"
"""

import sys

import matplotlib.pyplot as plt
import mne
import neurokit2 as nk


def open_one_participant(edf_path):
    # 1. Load the raw EDF file (MNE reads it without complaint, unlike pyedflib's strict parser)
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)

    sampling_rate = raw.info['sfreq']
    rsp_channel_index = raw.ch_names.index('Sensor-H:RSP')
    rsp_signal = raw.get_data(picks=[rsp_channel_index])[0]

    print(f'File: {edf_path}')
    print(f'Channels found: {raw.ch_names}')
    print(f'Sampling rate: {sampling_rate} Hz')
    print(f'Duration: {len(rsp_signal) / sampling_rate:.1f} s')
    print(f'Markers found: {list(raw.annotations.description)}')

    # 2. Visualize with NeuroKit2's own plotting function - raw signal only, nothing else.
    nk.signal_plot(rsp_signal, sampling_rate=sampling_rate)
    fig = plt.gcf()
    fig.set_size_inches(16, 4)
    return fig, rsp_signal, sampling_rate, raw


if __name__ == '__main__':
    edf_path = sys.argv[1]
    fig, rsp_signal, sampling_rate, raw = open_one_participant(edf_path)
    plt.show()  # opens an interactive window - zoom/pan with the mouse, close it to end the script
