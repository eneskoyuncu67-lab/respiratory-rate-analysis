# BART Respiratory Signal Pipeline

Recording and analysis pipeline for respiratory (chest-belt) signals collected during a BART (Balloon Analogue Risk Task) study, plus tools for reviewing and manually correcting automatically detected breaths.

This repository contains the **pipeline code only**. It does not include any participant data, recordings, or per-subject outputs — see [Setup](#setup) for how to point it at your own data.

## Contents

| File | Purpose |
| --- | --- |
| `config.py` | Central path configuration, read by every script below |
| `local_config.py.example` | Template for your own machine-specific paths (copy to `local_config.py`, which is gitignored) |
| `marker_parsing.py` | Shared logic for turning raw trigger markers into the 4 experimental block time ranges |
| `rsp_rate_by_block.py` | Main pipeline: computes automatic respiratory rate for every subject/block, writes the base CSV + summary plot + diagnostic plots |
| `manual_correct_single.py` | Interactive GUI to manually review/correct ONE chosen subject+block at a time |
| `manual_correct_peaks.py` | Interactive GUI that works through an automatic queue of every flagged (NA / low-confidence / fast) block |
| `rsp_summary_by_condition.py` | Aggregates per-block rates into a per-subject, per-condition (Certain/Uncertain) summary, including edit ratios |
| `rsp_correction_timeline.py` | Builds a clean, chronological CSV of every manual correction made |
| `watch_progress.py` | Quick progress/status check across all subjects and blocks |
| `rsp_full_recording_with_markers.py` | Visual QA: plots each subject's full raw recording with markers and resolved blocks shaded, for sanity-checking block detection |
| `combine_marker_plots_pdf.py` | Combines the QA plots above into a single multi-page PDF |
| `plot_with_markers.py` | Ad hoc: plot any single channel of one EDF file with trigger markers overlaid |
| `open_respiration.py` | Ad hoc: quickly view one subject's raw respiration trace, no processing |
| `rsp_breathmetrics_by_condition.m` / `plotViolin.m` / `rateFromSegments.m` | MATLAB alternative analysis using BreathMetrics instead of NeuroKit2 |

## Setup

1. Install Python dependencies:
   ```
   pip install -r requirements.txt
   ```
2. Either:
   - Drop your EDF recordings into `data/edf/` (created next to this README), and pipeline outputs will be written to `output/` — both work out of the box, or
   - Copy `local_config.py.example` to `local_config.py` (gitignored) and set `EDF_DIR` / `OUT_DIR` to your own locations, e.g. a separate data drive. Do the same with `local_config.m.example` → `local_config.m` for the MATLAB script.
3. Run the main pipeline:
   ```
   python rsp_rate_by_block.py
   ```
4. Manually review/correct flagged blocks (see [Manual correction](#manual-correction) below), then regenerate the summary:
   ```
   python rsp_summary_by_condition.py
   ```

### Third-party toolboxes

The MATLAB BreathMetrics pipeline expects [BreathMetrics](https://github.com/zelanolab/breathmetrics) to be available; point `BREATHMETRICS_ROOT` at it via `local_config.m` (see `local_config.m.example`).

---

## (a) Recording of respiratory signals

Respiration was recorded using a chest-expansion belt sampling the thoracic circumference at a native rate of 32 Hz, embedded as one channel within a multi-channel polysomnographic (EDF) recording alongside EEG and blood-volume-pulse channels sampled at higher native rates. Because the recording format stores all channels at a shared sampling rate, the respiration channel is internally upsampled to 256 Hz (matching the fastest channel in the file) when read via MNE-Python (Gramfort et al., 2013). This upsampling was verified to introduce no timing distortion: the true native 32 Hz samples, extracted directly from the raw EDF header and data records, were found to be identical to every 8th sample of the upsampled 256 Hz array (maximum absolute deviation = 0, correlation = 1.000).

Block onset and offset were marked by digital triggers logged at the start and end of each of the four experimental blocks (2 conditions × 2 repetitions), encoded as `block,condition` marker pairs (e.g. `1,2` = block 1, uncertain condition).

## (b) Analysis of respiratory data

Within each block, the respiratory signal was band-pass filtered (0.201–0.334 Hz, second-order Butterworth, following the filter architecture of the biosppy cleaning method implemented in NeuroKit2 — Makowski et al., 2021; Carreiras et al., 2015) and constant-detrended. The passband was derived from Rückert-Eheberg et al. (2025), who reported the 5th and 95th percentile respiratory rate (12.06 and 20.06 breaths/min, respectively) in a general-population sample of 2,224 adults (KORA-FF4 study), measured from resting electrocardiography.

This population-derived passband replaced an earlier, wider filter (0.05–3 Hz, NeuroKit2's default) after diagnostic inspection revealed that the wider default allowed non-respiratory signal components to be misclassified as breaths in a subset of recordings; this was confirmed via Welch's method (Welch, 1967), which showed affected blocks had peak spectral power outside the plausible respiratory range while unaffected blocks did not. Restricting the filter to the population-derived range eliminated the artifact in affected blocks while leaving unaffected blocks' estimates essentially unchanged.

Individual breaths were identified as local maxima in the filtered signal using a zero-crossing detection algorithm (Khodadad et al., 2018), excluding any candidate breath implying an instantaneous rate above 35 breaths per minute. Respiratory rate per block was calculated as 60 divided by the mean inter-peak interval in seconds.

Blocks were flagged for manual review if shorter than 60 seconds or containing fewer than 5 detected breaths, since such blocks provide insufficient data for a stable rate estimate; flagged blocks were retained (never excluded) but visually inspected against the underlying signal. Blocks exceeding 18 breaths/min were also included in review, following visual inspection of peak-interval distributions to distinguish genuine fast breathing from detection artifacts.

During review, individual detected peaks could be added or removed via a custom interactive correction interface (`manual_correct_single.py` / `manual_correct_peaks.py`) referenced against the filtered trace. For each manually reviewed block, an edit ratio was computed as the number of peaks added or removed (matched against the original automatic detection within a 0.5 s tolerance) divided by the number of peaks in that original detection, providing a quantitative index of how much manual intervention each estimate required. Both the original automatic detection and every manually corrected peak set are retained, so the full correction history of any block is reproducible.

## Manual correction

Two interactive tools are available, depending on your workflow:

- **`python manual_correct_single.py`** — prompts you for a subject and a block code (e.g. `08` then `1.1`), then opens a plot window for that one block only. Use this to jump straight to a specific subject/block.
- **`python manual_correct_peaks.py`** — automatically works through every block flagged NA / low-confidence / >18 bpm, one after another, in a queue.

Both tools:
- Print the last 5 correction actions from `correction_log.csv` on startup, so you always see where you left off.
- Save progress incrementally — closing partway through does not lose already-corrected blocks.
- Should **not** be run at the same time as each other (no file locking — running both simultaneously can silently overwrite each other's saves).

In the correction window:
- **Left-click near an existing detected peak** → removes that peak.
- **Left-click anywhere else on the trace** → adds a peak at the nearest local maximum within ±1s of the click.
- **Save & Next/Close** → saves your corrected peak list for this block.
- **Skip (keep original)** → leaves the automatic detection unchanged.

Block code convention: first digit = block number (1 or 2), second digit = condition (1 = certain, 2 = uncertain), e.g. `1.2` = 1st Uncertain.

## Outputs

All written to `OUT_DIR` (see [Setup](#setup)):

| File | Contents |
| --- | --- |
| `rsp_rate_by_block.csv` | Automatic (uncorrected) rates, all subjects/blocks |
| `rsp_rate_by_block_corrected.csv` | Same, updated live as you manually correct blocks |
| `rsp_summary_by_condition.csv` | Certain/Uncertain per-subject summary + edit ratios |
| `rsp_rate_by_block.png` | 4-column violin/dots summary plot |
| `manual_peaks_store.json` | Exact corrected peak timestamps per block (source of truth) |
| `correction_log.csv` | Chronological history of every save/discard action |
| `diagnostic_plots/` | One PNG per subject, all 4 blocks with detected peaks marked |

## Acknowledgments

Portions of this pipeline (refactoring for portability, documentation, and parts of the analysis code) were developed with AI assistance from Claude (Anthropic).

> Anthropic. (2026). *Claude (Sonnet 5)* [Large language model]. https://www.anthropic.com/claude

## References

- Carreiras, C., Alves, A. P., Lourenco, A., Canento, F., Silva, H., Fred, A., et al. (2015). BioSPPy - Biosignal Processing in Python. https://github.com/PIA-Group/BioSPPy/
- Gramfort, A., Luessi, M., Larson, E., Engemann, D. A., Strohmeier, D., Brodbeck, C., et al. (2013). MEG and EEG data analysis with MNE-Python. *Frontiers in Neuroscience*, 7, 267.
- Khodadad, D., Nordebo, S., Muller, B., Waldmann, A., Yerworth, R., Becher, T., et al. (2018). Optimized breath detection algorithm in electrical impedance tomography. *Physiological Measurement*, 39(9), 094001.
- Rückert-Eheberg, I.-M., Steger, A., Muller, A., Linkohr, B., et al. (2025). Respiratory rate and its associations with disease and lifestyle factors in the general population - results from the KORA-FF4 study. *PLOS ONE*, 20(3), e0318502.
- Makowski, D., Pham, T., Lau, Z. J., Brammer, J. C., Lespinasse, F., Pham, H., et al. (2021). NeuroKit2: A Python toolbox for neurophysiological signal processing. *Behavior Research Methods*, 53(4), 1689-1696.
- Welch, P. D. (1967). The use of fast Fourier transform for the estimation of power spectra: A method based on time averaging over short, modified periodograms. *IEEE Transactions on Audio and Electroacoustics*, 15(2), 70-73.
