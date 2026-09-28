"""
Per-participant respiratory rate using NeuroKit2, broken down per BLOCK,
displayed in the order: 1st Uncertain, 1st Certain, 2nd Uncertain,
2nd Certain.

NOTHING IS EXCLUDED. Every block that has a marker pair gets a rate
computed and plotted, including short/low-confidence ones (e.g. a block
whose start/end markers are only a few seconds apart, which is almost
certainly a double-tapped marker rather than a real short block - we
found one for zeynepacar's "1st Certain" block: 7.2s, 0 detected breaths).
Those are marked in RED instead of being dropped, so nothing is hidden -
you decide what to do with them, we don't decide for you.

A block is flagged low-confidence (red) if either:
  - duration < 60s (too little data for a stable average rate), or
  - fewer than 5 detected breaths (same reason)
These thresholds are just a display flag, not a filter - the numeric
value is still computed and still shown.

Cleaning uses a 2nd-order Butterworth bandpass filter (same architecture
as NeuroKit2's 'biosppy' cleaning method - Carreiras et al., 2015) plus
constant detrending, but with EVIDENCE-BASED cutoffs rather than
BioSPPy's software-default cutoffs (which were 0.1-0.35 Hz, i.e. 6-21
bpm - 6 bpm is bradypnea, not a real "normal adult" bound, and was never
actually backed by a citation, just inherited as a toolbox default).
The cutoffs used here instead come directly from a real measured
population: Ruckert-Eheberg et al. (2025, PLOS ONE) derived respiratory
rate from resting ECG in 2,224 general-population adults (KORA-FF4
study) and reported the 5th percentile at 12.06 breaths/min and the
95th percentile at 20.06 breaths/min. Converted to Hz (bpm/60), that's
RSP_LOWCUT_HZ = 0.2010 Hz and RSP_HIGHCUT_HZ = 0.3343 Hz below.

Peak detection still uses nk.rsp_peaks(method='biosppy') - that part of
the biosppy method is about extrema-finding logic (not filter cutoffs),
and was unaffected by the cleaning-cutoff change.

We originally reached the old 0.1-0.35 Hz passband by trial: NeuroKit2's
default cleaning method ('khodadad2018') is a very permissive 0.05-3 Hz
bandpass filter (deliberately covers 3-180 bpm), so ~0.5-1 Hz
contamination (cardiac bleed-through / movement artifact picked up by
the chestbelt) passed straight through untouched and got counted as
extra breaths. We confirmed this via Welch PSD: the "overcounted" blocks
(e.g. one subject's 1st Certain at 58.4 bpm, another's blocks at 35-42 bpm)
had their dominant spectral power sitting at ~0.5-1 Hz, well outside any
plausible respiratory band. Switching to a tight bandpass suppressed
that contamination at the source (58.4->16.2 bpm and
40.6->19.6 bpm at the time, respectively), which is what motivated tightening the
passband in the first place - the KORA-FF4-derived cutoffs used now are
a more defensible, citable version of that same fix.

We also checked whether the remaining >20 bpm blocks were still an
artifact (e.g. spurious doubling of single breaths into two peaks) by
looking at their peak-to-peak interval distributions. A doubling artifact
would show a bimodal distribution (a cluster of very short intervals mixed
with normal ones). Instead they were unimodal and continuous, with almost
no intervals below 2 seconds - i.e. real, not doubled - so those are left
as computed, not artificially suppressed.

Marker code convention (block-first, confirmed earlier in this project):
  first digit = block (1 or 2), second digit = condition (1=certain, 2=uncertain)
  each block is marked by the SAME code twice: once at start, once at end.
"""

# ---------------------------------------------------------------------------
# IMPORTS - external libraries this script depends on.
# ---------------------------------------------------------------------------
import glob   # glob.glob() finds files matching a wildcard pattern, e.g. "*.edf"
import os     # for building file paths (os.path.join) and creating folders
import re     # regular expressions - used to recognize marker text like "1,2"

import matplotlib
matplotlib.use('Agg')  # tells matplotlib to render to image files, NOT open a
                        # window - required because this script runs headless
                        # (no display) when generating plots in bulk. The
                        # interactive correction tools use 'TkAgg' instead,
                        # because THEY need a real window you can click on.
import matplotlib.pyplot as plt  # the actual plotting functions (plot, scatter, etc.)
import mne         # reads the EDF recordings (signal channels + trigger markers)
import neurokit2 as nk  # the physiological-signal library: still used for
                         # nk.rsp_peaks(method='biosppy') peak detection
from neurokit2.signal import signal_detrend  # same detrending NeuroKit2's
                         # own biosppy cleaning method uses internally
import numpy as np      # array math: means, differences between numbers, etc.
import pandas as pd     # builds the final results table and writes it to CSV
from scipy.signal import butter, filtfilt  # build our own bandpass filter,
                         # since NeuroKit2's biosppy method hardcodes its
                         # cutoffs and doesn't let us pass different ones


# ---------------------------------------------------------------------------
# CONFIGURATION - paths and constants used throughout the script.
# ---------------------------------------------------------------------------

# Where the raw recordings live and where all outputs (CSV table, summary
# plot, per-subject diagnostic plots) get written - see config.py to point
# these at your own data location.
from config import EDF_DIR, OUT_DIR

# Subfolder specifically for the per-subject diagnostic plots (one PNG per
# subject, showing all 4 of their blocks side by side).
DIAG_DIR = os.path.join(OUT_DIR, 'diagnostic_plots')

# Create these folders if they don't already exist. exist_ok=True means "don't
# raise an error if the folder is already there" - so this is safe to re-run.
os.makedirs(DIAG_DIR, exist_ok=True)

# A "regular expression" pattern that matches marker text shaped like "1,2" or
# "2,1" - i.e. exactly one digit (1 or 2), a comma, then exactly one more
# digit (1 or 2). Anything that doesn't look like this (noise markers,
# stray text) will simply not match and gets skipped later.
CODE_RE = re.compile(r'^([12]),([12])$')

# The 4 blocks, in the exact left-to-right order we want them to appear in
# the final plot: 1st Uncertain, 1st Certain, 2nd Uncertain, 2nd Certain.
# Each tuple is (block_digit, condition_digit, human-readable label).
# Marker convention: first digit = block number (1 or 2), second digit =
# condition (1 = certain, 2 = uncertain). So ('1','2', ...) means "block 1,
# condition 2 (uncertain)".
BLOCK_KEYS = [('1', '2', '1st Uncertain'), ('1', '1', '1st Certain'),
              ('2', '2', '2nd Uncertain'), ('2', '1', '2nd Certain')]

# A block gets flagged "low confidence" (shown in red) if it's shorter than
# this many seconds - too little data to trust an average rate from.
MIN_DURATION_S = 60

# ...or if fewer than this many breaths were detected in it - same reason,
# just measured in breath-count instead of time.
MIN_PEAKS = 5

# Bandpass filter cutoffs, in Hz, converted from breaths/min by dividing by
# 60 (60 seconds in a minute). Source: Ruckert-Eheberg et al. (2025, PLOS
# ONE), 5th and 95th percentile respiratory rate in a 2,224-person general-
# population sample (KORA-FF4 study), derived from resting ECG.
RSP_LOWCUT_HZ = 12.06 / 60   # = 0.2010 Hz (5th percentile, 12.06 breaths/min)
RSP_HIGHCUT_HZ = 20.06 / 60  # = 0.3343 Hz (95th percentile, 20.06 breaths/min)
FILTER_ORDER = 2  # matches BioSPPy's own filter order (Carreiras et al., 2015)


def clean_rsp(seg, sfreq):
    """
    Band-pass filters + detrends one block's raw RSP segment, using the
    same filter ARCHITECTURE as NeuroKit2's biosppy cleaning method
    (2nd-order Butterworth bandpass, then constant detrending) but with
    our own evidence-based cutoffs (RSP_LOWCUT_HZ/RSP_HIGHCUT_HZ above)
    instead of BioSPPy's hardcoded 0.1-0.35 Hz defaults, which we can't
    change by passing a parameter - NeuroKit2's biosppy method has those
    numbers written directly into its source code.
    """
    # scipy.signal.butter wants cutoff frequencies normalized to the
    # Nyquist frequency (half the sampling rate), not raw Hz - this line
    # does that conversion for both the low and high cutoff at once.
    freq = 2 * np.array([RSP_LOWCUT_HZ, RSP_HIGHCUT_HZ]) / sfreq

    # Design the filter (returns its coefficients, b and a) ...
    b, a = butter(N=FILTER_ORDER, Wn=freq, btype='bandpass')
    # ...then apply it. filtfilt runs the filter forward AND backward,
    # which cancels out the phase-shift/lag a normal one-directional
    # filter would introduce - so breath timing stays accurate.
    filtered = filtfilt(b, a, seg)

    # Remove any constant (DC) offset in the signal - same final step
    # NeuroKit2's own biosppy cleaning method performs.
    return signal_detrend(filtered, order=0)


# ---------------------------------------------------------------------------
# Block-boundary detection now lives in marker_parsing.py, shared by all
# three scripts (this one + both correction tools) so they can never
# disagree with each other. It combines a general parser (accepts comma
# OR period marker codes, debounces double-taps, falls back to a
# single-occurrence convention where needed) with a small table of
# hand-diagnosed per-subject overrides for markers that were recorded
# with the wrong code or are missing entirely - see that file for the
# full reasoning behind each override.
# ---------------------------------------------------------------------------
from marker_parsing import get_block_segments


# ---------------------------------------------------------------------------
# FUNCTION: given one block's raw signal chunk, clean it and detect
# individual breaths in it, then compute a rate.
# ---------------------------------------------------------------------------
def peaks_for_segment(sig, sfreq, start_s, end_s):
    """
    Inputs:
      sig      - the FULL respiration signal for this participant (one long
                 array covering the entire recording)
      sfreq    - sampling rate in Hz (samples per second)
      start_s, end_s - the block's start/end time in seconds, from
                 get_block_segments()

    Returns a 4-tuple (rate_bpm, peaks, cleaned, low_confidence):
      rate_bpm       - breaths per minute, or NaN ("not a number") if fewer
                        than 2 breaths were found (you need at least 2
                        breaths to measure the time between them)
      peaks           - the sample-index position of every detected breath,
                        counted from the start of THIS block (not the whole
                        recording) - used later to draw red dots on the plot
      cleaned         - the cleaned (filtered) version of this block's
                        signal, also used for plotting
      low_confidence  - True if this block is short or has few breaths,
                        meaning the rate above should be treated with
                        caution - but it's still returned, never hidden
    """
    # Convert the block's start/end times (in seconds) into sample-index
    # positions within the full signal array, then slice out just this
    # block's chunk of data.
    i0, i1 = int(start_s * sfreq), int(end_s * sfreq)
    seg = sig[i0:i1]
    duration_s = end_s - start_s

    # If there are fewer than 10 samples total, this block is essentially
    # empty (almost certainly a double-tapped marker, not a real block) -
    # there's nothing meaningful to filter or detect peaks in, so bail out
    # immediately with NaN and flag it low-confidence.
    if len(seg) < 10:
        return np.nan, np.array([], dtype=int), seg, True

    try:
        # Step 1: CLEAN the signal - removes noise/drift and band-passes it
        # to RSP_LOWCUT_HZ-RSP_HIGHCUT_HZ (evidence-based cutoffs, see the
        # big comment at the top of this file for the citation and why).
        cleaned = clean_rsp(seg, sfreq)

        # Step 2: DETECT PEAKS - finds each individual breath (specifically,
        # the moment of peak inhalation) in the cleaned signal. The biosppy
        # method also automatically rejects any pair of peaks that are so
        # close together they'd imply an impossible rate above 35 bpm.
        _, info = nk.rsp_peaks(cleaned, sampling_rate=sfreq, method='biosppy')
        peaks = info['RSP_Peaks']  # array of sample-index positions of each breath peak
    except Exception:
        # If anything goes wrong during cleaning/detection (e.g. a weird
        # edge case NeuroKit2 can't handle), treat it the same as "nothing
        # found" rather than crashing the whole script.
        return np.nan, np.array([], dtype=int), seg, True

    # Decide whether to flag this block as low-confidence: either it's too
    # short in time, or too few breaths were found in it.
    low_confidence = (duration_s < MIN_DURATION_S) or (len(peaks) < MIN_PEAKS)

    if len(peaks) < 2:
        # You need at least 2 breath peaks to measure the time GAP between
        # them - with 0 or 1 peaks there's no interval to measure, so no
        # rate can be computed at all. This isn't a choice, it's just math.
        return np.nan, peaks, cleaned, True

    # Rate calculation:
    #   np.diff(peaks)       -> the number of SAMPLES between each
    #                            consecutive pair of breath peaks
    #   np.diff(peaks)/sfreq -> convert that from samples to SECONDS
    #   np.mean(...)         -> average seconds-per-breath across the block
    #   60 / average         -> convert "seconds per breath" into
    #                            "breaths per minute" (60 seconds in a minute)
    rate = 60 / np.mean(np.diff(peaks) / sfreq)
    return rate, peaks, cleaned, low_confidence


# ---------------------------------------------------------------------------
# MAIN LOOP - go through every recording, extract all 4 blocks from each,
# and build up the results table row by row.
# ---------------------------------------------------------------------------

# Find every .edf file in the folder and sort them alphabetically (which,
# since filenames start with the 2-digit participant number, also puts them
# in the correct participant order).
files = sorted(glob.glob(os.path.join(EDF_DIR, '*.edf')))

rows = []  # will collect one dictionary per participant, becomes the final table

for f in files:
    name = os.path.basename(f)  # just "05_subjectname.edf", not the full folder path
    try:
        # Load the whole recording into memory. preload=True means "read
        # the actual signal data now" (not just the header info).
        # verbose=False silences MNE's routine loading messages.
        raw = mne.io.read_raw_edf(f, preload=True, verbose=False)
        sfreq = raw.info['sfreq']  # sampling rate in Hz for this recording

        # Pull out just the respiration channel's data as a plain 1-D array.
        # raw.get_data(picks=[...]) returns a 2-D array (channels x samples)
        # even when you only ask for one channel, so [0] grabs that single row.
        sig = raw.get_data(picks=[raw.ch_names.index('Sensor-H:RSP')])[0]

        # Find this participant's 4 block time-windows from their markers.
        segments = get_block_segments(raw, name)

        row = {'file': name}  # this participant's row in the final table, so far just their filename
        diag_segments = []  # will hold data needed to draw this participant's diagnostic plot

        # Go through the 4 blocks in the display order we defined earlier.
        for block, cond, label in BLOCK_KEYS:
            key = f'{block}_{cond}'
            if key not in segments:
                # This participant has no marker pair for this particular
                # block (e.g. their trigger scheme was unusable) - record it
                # as missing data and move on to the next block.
                row[label] = np.nan
                row[f'{label}_low_conf'] = False
                continue

            start_s, end_s = segments[key]
            # Do the actual cleaning + peak detection + rate calculation
            # for this one block.
            rate, peaks, cleaned, low_conf = peaks_for_segment(sig, sfreq, start_s, end_s)

            row[label] = rate
            row[f'{label}_low_conf'] = low_conf
            # Remember everything needed to plot this block later.
            diag_segments.append((label, start_s, end_s, peaks, cleaned, low_conf))

        rows.append(row)  # add this participant's finished row to the table

        # Build a one-line human-readable summary of this participant's 4
        # rates to print to the console, e.g. "1st Uncertain=17.7  1st
        # Certain=19.6(LOW-CONF)  ...". The bit inside the f-string picks
        # between showing the number or "NA" depending on whether it's NaN.
        vals = '  '.join(
            f'{label}={row[label]:.1f}{"(LOW-CONF)" if row[f"{label}_low_conf"] else ""}'
            if not np.isnan(row[label]) else f'{label}=NA'
            for _, _, label in BLOCK_KEYS
        )
        print(f'{name:30s} {vals}')

        # -----------------------------------------------------------------
        # Per-participant diagnostic plot: one row of panels, one panel per
        # block that had data, showing the cleaned signal with every
        # detected breath marked as a red dot. Low-confidence blocks get a
        # red title so they stand out visually when you're scanning through
        # many of these PNGs.
        # -----------------------------------------------------------------
        if diag_segments:  # only bother plotting if at least one block had data
            # Create one row of side-by-side subplots, one per block found.
            # figsize scales with how many panels there are so they don't
            # get squeezed. squeeze=False keeps 'axes' as a consistent 2-D
            # array even if there's only 1 panel, so the indexing below
            # always works the same way.
            fig, axes = plt.subplots(1, len(diag_segments), figsize=(5 * len(diag_segments), 4),
                                      squeeze=False)

            # zip() pairs up each subplot axis with its corresponding
            # block's data, so we can draw each block into its own panel.
            for ax, (label, start_s, end_s, peaks, cleaned, low_conf) in zip(axes[0], diag_segments):
                t = np.arange(len(cleaned)) / sfreq  # time axis in seconds, starting at 0 for this block
                ax.plot(t, cleaned, 'k-', linewidth=0.6)  # the cleaned signal itself, thin black line
                ax.scatter(t[peaks], cleaned[peaks], color='red', s=20, zorder=3)  # red dots on each detected breath
                title_color = 'red' if low_conf else 'black'
                ax.set_title(f'{label}{" [LOW CONF]" if low_conf else ""}\n'
                              f'({len(peaks)} breaths, {(end_s-start_s)/60:.1f} min)',
                              color=title_color)
                ax.set_xlabel('Time in block (s)')

            axes[0][0].set_ylabel('Cleaned RSP')  # only the leftmost panel needs a y-axis label
            fig.suptitle(name, fontsize=10)  # participant's filename as the overall figure title
            plt.tight_layout()  # auto-adjusts spacing so labels/titles don't overlap
            # Save as "<participant name without .edf>.png" into the diagnostic-plots folder.
            plt.savefig(os.path.join(DIAG_DIR, f'{os.path.splitext(name)[0]}.png'), dpi=120)
            plt.close(fig)  # free up memory - important when looping over 50 participants

    except Exception as e:
        # If anything about this participant's file fails entirely (e.g.
        # corrupted file, unexpected format), don't let it crash the whole
        # batch - just report it and move on to the next participant.
        print(f'{name}: FAILED {e}')

# ---------------------------------------------------------------------------
# Turn the list of per-participant dictionaries into a proper table and
# save it as a CSV file you can open in Excel or load back into Python.
# ---------------------------------------------------------------------------
df = pd.DataFrame(rows)
csv_path = os.path.join(OUT_DIR, 'rsp_rate_by_block.csv')
df.to_csv(csv_path, index=False)
print(f'\nSaved: {csv_path}')
print(f'Saved diagnostic plots to: {DIAG_DIR}')

# ---------------------------------------------------------------------------
# FINAL SUMMARY PLOT: one violin (distribution shape) + a swarm of dots per
# block, in the order 1st Uncertain / 1st Certain / 2nd Uncertain / 2nd
# Certain, with faint lines connecting each participant's 4 points. Every
# value is plotted - nothing excluded - but low-confidence points are red.
# ---------------------------------------------------------------------------
labels = [label for _, _, label in BLOCK_KEYS]  # just the 4 column names, in order
fig, ax = plt.subplots(figsize=(10, 8))

# Gather each column's non-missing values into a list of arrays - one array
# per block - to feed into the violin plot.
data = [df[label].dropna().values for label in labels]
positions = list(range(1, len(labels) + 1))  # x-axis positions 1,2,3,4 for the 4 blocks

# Draw the violin shapes. If a column has fewer than 2 real values,
# matplotlib can't compute a density shape for it, so we substitute a
# [NaN, NaN] placeholder for that column to avoid crashing.
parts = ax.violinplot([d if len(d) >= 2 else [np.nan, np.nan] for d in data],
                       positions=positions, showmeans=True, showextrema=True)
for pc in parts['bodies']:
    pc.set_facecolor('lightgray')
    pc.set_alpha(0.5)

# Scatter each participant's individual value as a dot on top of the violin,
# with small random left-right jitter so overlapping dots are still visible.
rng = np.random.default_rng(0)  # fixed seed -> same jitter pattern every time you rerun this
for pos, label in zip(positions, labels):
    # Grab this column's values AND its matching low-confidence flags,
    # dropping rows where the value itself is missing.
    sub = df[[label, f'{label}_low_conf']].dropna(subset=[label])
    jitter = rng.uniform(-0.05, 0.05, size=len(sub))
    # Pick red for low-confidence points, black otherwise - one color per dot.
    colors = np.where(sub[f'{label}_low_conf'].values, 'red', 'black')
    ax.scatter(np.full(len(sub), pos) + jitter, sub[label].values,
               alpha=0.7, s=22, c=colors, zorder=3)

# Draw a faint gray line connecting each participant's 4 points, but ONLY
# for participants who have a real value in all 4 columns (dropna with
# subset=labels keeps only fully-complete rows).
complete = df.dropna(subset=labels)
for _, r in complete.iterrows():
    ax.plot(positions, [r[label] for label in labels], color='gray', alpha=0.15, linewidth=0.8, zorder=1)

# Shade the 10-18 bpm range in light green as a rough "plausible adult
# resting rate" reference band, purely visual - doesn't affect any numbers.
ax.axhspan(10, 18, color='green', alpha=0.08, zorder=0)

# Add two invisible scatter points purely to generate legend entries
# explaining what black vs red dots mean.
ax.scatter([], [], c='black', s=22, label='rate (\u226560s block, \u22655 breaths)')
ax.scatter([], [], c='red', s=22, label=f'low confidence (<{MIN_DURATION_S}s block or <{MIN_PEAKS} breaths)')
ax.legend(loc='upper right', fontsize=8)

ax.set_xticks(positions)
ax.set_xticklabels(labels)
ax.set_ylabel('Respiratory rate (breaths/min), NeuroKit2')
ax.set_title(f'Respiratory rate per block (n={len(df)}), nothing excluded')
plt.tight_layout()
out_plot = os.path.join(OUT_DIR, 'rsp_rate_by_block.png')
plt.savefig(out_plot, dpi=140)
print(f'Saved: {out_plot}')
