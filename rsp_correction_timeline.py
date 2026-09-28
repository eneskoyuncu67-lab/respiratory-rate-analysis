"""
Builds a clean timeline CSV of every block you've actually corrected (not
discards, not internal bookkeeping columns) - one row per correction, in
the order you did them, with: timestamp, subject, block, resulting
respiratory rate, and the edit ratio for that correction.

Reads from correction_log.csv (the raw event log both correction tools
write to) and filters/reshapes it into this cleaner view.
"""
import os

import pandas as pd

from config import OUT_DIR
LOG_PATH = os.path.join(OUT_DIR, 'correction_log.csv')
OUT_PATH = os.path.join(OUT_DIR, 'rsp_correction_timeline.csv')

if not os.path.exists(LOG_PATH):
    print(f'No correction history yet at {LOG_PATH} - correct at least one block first.')
    raise SystemExit

log = pd.read_csv(LOG_PATH)

# Only keep actual saves - discards have no rate/ratio and aren't a
# "correction", just a record that you looked at something and left it.
saved = log[log['action'] == 'saved'].copy()

if saved.empty:
    print('No saved corrections yet (only discards, or nothing at all).')
    raise SystemExit

# Keep just the columns relevant to "rate + correction ratio per block",
# renamed to be self-explanatory as a standalone table.
timeline = saved[['timestamp', 'file', 'label', 'rate_bpm', 'n_original',
                   'n_added', 'n_removed', 'edit_ratio']].copy()
timeline.columns = ['timestamp', 'subject', 'block', 'respiratory_rate_bpm',
                     'n_original_peaks', 'n_added', 'n_removed', 'edit_ratio']

# Sort chronologically (oldest first) so it reads top-to-bottom as a
# timeline of your actual work session order - if you corrected the same
# block twice (discarded then redone, or revised later), both entries
# stay, in order, so you can see the history of that block too.
timeline = timeline.sort_values('timestamp').reset_index(drop=True)

timeline.to_csv(OUT_PATH, index=False)

print(f'Saved: {OUT_PATH}')
print(f'{len(timeline)} correction(s) recorded.\n')
print(timeline.to_string(index=False))
