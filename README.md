# BART Respiratory Signal Pipeline

This pipeline turns chest-belt respiration recordings into a respiratory rate for each experimental block. It was built for a Balloon Analogue Risk Task (BART) design. It also includes interactive tools to review detected breaths and correct them by hand.

The repository holds the pipeline code only. It contains no participant data, recordings or per-subject outputs. See [Setup](#setup) to point it at your own data.

## Contents

| File | Purpose |
| --- | --- |
| `config.py` | Central path settings used by every script |
| `local_config.py.example` | Template for machine-specific paths. Copy it to `local_config.py`, which git ignores |
| `marker_parsing.py` | Turns raw trigger markers into the time ranges of the four experimental blocks |
| `rsp_rate_by_block.py` | Main pipeline. Computes the automatic respiratory rate for every subject and block. Writes the base CSV, a summary plot and diagnostic plots |
| `manual_correct_single.py` | Interactive tool to review and correct one chosen subject and block |
| `manual_correct_peaks.py` | Interactive tool that steps through a queue of all flagged blocks (NA, low confidence or fast) |
| `rsp_summary_by_condition.py` | Summarises block rates per subject and condition (Certain or Uncertain). Includes edit ratios |
| `rsp_correction_timeline.py` | Builds a chronological CSV of every manual correction |
| `watch_progress.py` | Shows correction progress across all subjects and blocks |
| `rsp_full_recording_with_markers.py` | Quality check. Plots each full raw recording with its markers and shades the detected blocks |
| `combine_marker_plots_pdf.py` | Combines these quality-check plots into one multi-page PDF |
| `plot_with_markers.py` | Plots any single channel of one EDF file with its trigger markers |
| `open_respiration.py` | Shows one subject's raw respiration trace without processing |
| `rsp_breathmetrics_by_condition.m` / `plotViolin.m` / `rateFromSegments.m` | Alternative MATLAB analysis that uses BreathMetrics instead of NeuroKit2 |

## Setup

1. Install the Python dependencies:
   ```
   pip install -r requirements.txt
   ```
2. Choose where data and outputs live. There are two options:
   - Put your EDF recordings in `data/edf/` next to this README. Outputs then go to `output/`. This works without further setup.
   - Copy `local_config.py.example` to `local_config.py` and set `EDF_DIR` and `OUT_DIR` to your own folders. For the MATLAB script do the same with `local_config.m.example` and `local_config.m`.
3. Run the main pipeline:
   ```
   python rsp_rate_by_block.py
   ```
4. Review and correct the flagged blocks (see [Manual correction](#manual-correction)). Then rebuild the summary:
   ```
   python rsp_summary_by_condition.py
   ```

### Third-party toolboxes

The MATLAB analysis needs [BreathMetrics](https://github.com/zelanolab/breathmetrics). Set `BREATHMETRICS_ROOT` to its location in `local_config.m` (see `local_config.m.example`).

---

## How the pipeline works

### Input signals

The pipeline expects respiration from a chest-expansion belt. The belt signal is one channel of a multi-channel EDF recording. It was developed for a belt sampled at 32 Hz next to EEG and blood-volume-pulse channels with higher sampling rates.

The recording format stores all channels at one shared sampling rate. When MNE-Python (Gramfort et al., 2013) reads the file the respiration channel is therefore upsampled to 256 Hz. This matches the fastest channel in the file. The upsampling does not distort timing. The native 32 Hz samples read directly from the raw EDF records are identical to every 8th sample of the 256 Hz array (maximum absolute deviation 0, correlation 1.000).

Digital triggers mark the start and end of each block. The task has four blocks made of two conditions with two repetitions each. Markers come as `block,condition` pairs. For example `1,2` means block 1 in the uncertain condition.

### Filtering

Within each block the signal is band-pass filtered from 0.201 to 0.334 Hz and then constant-detrended. The filter is a second-order Butterworth. It follows the filter design of the biosppy cleaning method in NeuroKit2 (Makowski et al., 2021; Carreiras et al., 2015).

The passband comes from Rückert-Eheberg et al. (2025). They report the 5th and 95th percentiles of respiratory rate as 12.06 and 20.06 breaths per minute. These values were measured from resting electrocardiography in 2,224 adults of the general-population KORA-FF4 study.

The pipeline does not use the NeuroKit2 default passband of 0.05 to 3 Hz. This wider band can let non-respiratory components through as breaths. Welch's method (Welch, 1967) makes this visible. Affected blocks show their spectral peak outside the plausible respiratory range. In testing the narrower band removed these artifacts. Estimates for unaffected blocks stayed essentially the same.

### Breath detection

Breaths are local maxima in the filtered signal. A zero-crossing detection algorithm finds them (Khodadad et al., 2018). Any candidate that implies an instantaneous rate above 35 breaths per minute is discarded. The respiratory rate of a block is 60 divided by the mean inter-peak interval in seconds.

### Flagging blocks for review

A block is flagged for manual review if it is shorter than 60 s or has fewer than 5 detected breaths. Such blocks hold too little data for a stable rate. Blocks faster than 18 breaths per minute are flagged as well. This threshold was set from the peak-interval distributions. It helps separate genuine fast breathing from detection artifacts. Flagged blocks are never excluded. They are checked by eye against the underlying signal.

### Correction and edit ratio

The correction tools let you add or remove single peaks on top of the filtered trace. Each reviewed block gets an edit ratio. It is the number of peaks added or removed divided by the number of automatically detected peaks. Corrected peaks are matched to the automatic ones within a 0.5 s tolerance. The ratio shows how much manual work each estimate needed.

The pipeline keeps both the automatic detection and every corrected peak set. The full correction history of any block can therefore be reproduced.

## Manual correction

Two interactive tools cover different workflows:

- **`python manual_correct_single.py`** asks for a subject and a block code (e.g. `08` then `1.1`). It then opens a plot of that one block. Use it to jump straight to a specific block.
- **`python manual_correct_peaks.py`** steps through every block flagged as NA, low confidence or above 18 bpm.

Both tools:
- Show the last 5 actions from `correction_log.csv` at startup. You always see where you left off.
- Save after each block. Closing early keeps the blocks you have already corrected.
- Must **not** run at the same time. There is no file locking. Two running tools can silently overwrite each other's saves.

In the correction window:
- **Left-click near a detected peak** to remove it.
- **Left-click elsewhere on the trace** to add a peak at the nearest local maximum within ±1 s.
- **Save & Next/Close** saves the corrected peaks for this block.
- **Skip (keep original)** keeps the automatic detection.

Block codes work like this. The first digit is the block number (1 or 2). The second digit is the condition (1 = certain, 2 = uncertain). For example `1.2` is the first uncertain block.

## Outputs

All files are written to `OUT_DIR` (see [Setup](#setup)).

| File | Contents |
| --- | --- |
| `rsp_rate_by_block.csv` | Automatic (uncorrected) rates for all subjects and blocks |
| `rsp_rate_by_block_corrected.csv` | The same table, updated as you correct blocks |
| `rsp_summary_by_condition.csv` | Per-subject summary by condition with edit ratios |
| `rsp_rate_by_block.png` | Summary plot with violins and dots for the four blocks |
| `manual_peaks_store.json` | Exact corrected peak times per block. This file is the source of truth |
| `correction_log.csv` | History of every save and discard |
| `diagnostic_plots/` | One PNG per subject with all four blocks and their detected peaks |

## Acknowledgments

Parts of this pipeline were developed with AI assistance from Claude (Anthropic). This covers refactoring for portability, documentation and parts of the analysis code.

> Anthropic. (2026). *Claude (Sonnet 5)* [Large language model]. https://www.anthropic.com/claude

## References

- Carreiras, C., Alves, A. P., Lourenco, A., Canento, F., Silva, H., Fred, A., et al. (2015). BioSPPy - Biosignal Processing in Python. https://github.com/PIA-Group/BioSPPy/
- Gramfort, A., Luessi, M., Larson, E., Engemann, D. A., Strohmeier, D., Brodbeck, C., et al. (2013). MEG and EEG data analysis with MNE-Python. *Frontiers in Neuroscience*, 7, 267. https://doi.org/10.3389/fnins.2013.00267
- Khodadad, D., Nordebo, S., Muller, B., Waldmann, A., Yerworth, R., Becher, T., et al. (2018). Optimized breath detection algorithm in electrical impedance tomography. *Physiological Measurement*, 39(9), 094001. https://doi.org/10.1088/1361-6579/aad7e6
- Rückert-Eheberg, I.-M., Steger, A., Muller, A., Linkohr, B., et al. (2025). Respiratory rate and its associations with disease and lifestyle factors in the general population - results from the KORA-FF4 study. *PLOS ONE*, 20(3), e0318502. https://doi.org/10.1371/journal.pone.0318502
- Makowski, D., Pham, T., Lau, Z. J., Brammer, J. C., Lespinasse, F., Pham, H., et al. (2021). NeuroKit2: A Python toolbox for neurophysiological signal processing. *Behavior Research Methods*, 53(4), 1689-1696. https://doi.org/10.3758/s13428-020-01516-y 
