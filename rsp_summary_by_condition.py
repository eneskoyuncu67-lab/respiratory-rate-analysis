"""
Collapses the 4 per-block results (1st Certain, 1st Uncertain, 2nd Certain,
2nd Uncertain) down to 2 per-CONDITION numbers per subject (Certain,
Uncertain), by averaging that condition's two blocks - and alongside each
rate, reports how much manual correction went into it.

Reads rsp_rate_by_block_corrected.csv if it exists (so any manual
corrections you've made are reflected), otherwise falls back to the
uncorrected rsp_rate_by_block.csv.

For each subject and condition, reports:
  {condition}_rate         - mean bpm across that condition's 2 blocks
                              (using whatever value is currently in the
                              CSV for each block - corrected if you've
                              corrected it, original NeuroKit2 output if not)
  {condition}_n_reviewed    - how many of the 2 contributing blocks have
                              actually been manually reviewed (0, 1, or 2)
  {condition}_edit_ratio    - mean edit_ratio across only the blocks that
                              HAVE been reviewed (ignores blocks not yet
                              touched, rather than treating them as 0 -
                              "not reviewed" and "reviewed, 0 edits" are
                              different things and shouldn't be conflated)
"""
import os

import numpy as np
import pandas as pd

from config import OUT_DIR
CORRECTED_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block_corrected.csv')
BASE_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block.csv')

# Prefer the corrected file (has your manual edits); fall back to the raw
# automatic output if no corrections have been made yet at all.
if os.path.exists(CORRECTED_CSV):
    df = pd.read_csv(CORRECTED_CSV)
    print(f'Loaded (with corrections): {CORRECTED_CSV}')
else:
    df = pd.read_csv(BASE_CSV)
    print(f'Loaded (no corrections yet): {BASE_CSV}')

# Which two block-columns feed into each condition.
CONDITIONS = {
    'Certain': ['1st Certain', '2nd Certain'],
    'Uncertain': ['1st Uncertain', '2nd Uncertain'],
}

summary_rows = []
for _, row in df.iterrows():
    out = {'file': row['file']}
    for cond_name, block_labels in CONDITIONS.items():
        rates = []
        edit_ratios_reviewed = []
        n_reviewed = 0
        for label in block_labels:
            val = row.get(label, np.nan)
            if pd.notna(val):
                rates.append(val)
            # '_corrected' column only exists if at least one block anywhere
            # in the whole dataset has been through the correction tool -
            # use .get(..., False) so this doesn't crash before that.
            was_reviewed = bool(row.get(f'{label}_corrected', False))
            if was_reviewed:
                n_reviewed += 1
                er = row.get(f'{label}_edit_ratio', np.nan)
                if pd.notna(er):
                    edit_ratios_reviewed.append(er)

        out[f'{cond_name}_rate'] = np.mean(rates) if rates else np.nan
        out[f'{cond_name}_n_reviewed'] = n_reviewed  # out of 2 blocks
        out[f'{cond_name}_edit_ratio'] = (np.mean(edit_ratios_reviewed)
                                           if edit_ratios_reviewed else np.nan)
    summary_rows.append(out)

summary = pd.DataFrame(summary_rows)
out_path = os.path.join(OUT_DIR, 'rsp_summary_by_condition.csv')
summary.to_csv(out_path, index=False)

print(f'Saved: {out_path}')
print()
print(summary.to_string(index=False))

# ---------------------------------------------------------------------------
# Overall correction-progress stats, printed for a quick sanity check.
# ---------------------------------------------------------------------------
block_labels = ['1st Certain', '1st Uncertain', '2nd Certain', '2nd Uncertain']
total_blocks = 0
reviewed_blocks = 0
all_edit_ratios = []
for label in block_labels:
    col_corrected = f'{label}_corrected'
    col_ratio = f'{label}_edit_ratio'
    if col_corrected not in df.columns:
        continue  # nobody has corrected any block of this type yet
    total_blocks += df[label].notna().sum()  # only count blocks that have a real value at all
    reviewed_mask = df[col_corrected] == True
    reviewed_blocks += reviewed_mask.sum()
    if col_ratio in df.columns:
        all_edit_ratios.extend(df.loc[reviewed_mask, col_ratio].dropna().tolist())

print()
print(f'Overall: {reviewed_blocks} / {total_blocks} existing blocks manually reviewed so far.')
if all_edit_ratios:
    print(f'Mean edit ratio across reviewed blocks: {np.mean(all_edit_ratios):.3f}')
else:
    print('No reviewed blocks with a recorded edit ratio yet.')
