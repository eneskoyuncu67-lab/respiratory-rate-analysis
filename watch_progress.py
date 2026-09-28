"""
Live progress dashboard - run this in a SEPARATE terminal window (e.g. on
your second screen) while you run manual_correct_single.py or
manual_correct_peaks.py in your main window. This script does not
interact with anything - it just re-reads the log/CSV files every few
seconds and reprints an up-to-date summary, so you can watch your
progress update in real time as you save corrections elsewhere.

Run this locally:
    D:\\bart_v1\\neurokit\\python.exe D:\\bart_v1\\watch_progress.py

Press Ctrl+C to stop watching.
"""
import os
import time

import numpy as np
import pandas as pd

from config import OUT_DIR
LOG_PATH = os.path.join(OUT_DIR, 'correction_log.csv')
CORRECTED_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block_corrected.csv')
BASE_CSV = os.path.join(OUT_DIR, 'rsp_rate_by_block.csv')

REFRESH_SECONDS = 5
LOG_ENTRIES_TO_SHOW = 15

BLOCK_LABELS = ['1st Uncertain', '1st Certain', '2nd Uncertain', '2nd Certain']


def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')


def render_once():
    clear_screen()
    print('=' * 78)
    print(f'  LIVE CORRECTION PROGRESS - refreshing every {REFRESH_SECONDS}s '
          f'(Ctrl+C to stop)')
    print(f'  {time.strftime("%Y-%m-%d %H:%M:%S")}')
    print('=' * 78)
    print()

    # -----------------------------------------------------------------
    # Recent history (timestamps) - the "timestamp folder" info: what
    # you've done most recently, in order.
    # -----------------------------------------------------------------
    if os.path.exists(LOG_PATH):
        log = pd.read_csv(LOG_PATH)
        print(f'Last {min(LOG_ENTRIES_TO_SHOW, len(log))} action(s):')
        for _, r in log.tail(LOG_ENTRIES_TO_SHOW).iterrows():
            rate = r.get('rate_bpm', '')
            rate_str = f'{rate:.1f}' if pd.notna(rate) and rate != '' else '-'
            print(f'  {r["timestamp"]}  [{r["tool"]:6s}]  {r["file"]:28s} '
                  f'{r["label"]:15s} {r["action"]:9s} rate={rate_str}')
        print()
    else:
        print('No correction history yet.\n')

    # -----------------------------------------------------------------
    # Overall progress: how many blocks have been reviewed so far.
    # -----------------------------------------------------------------
    csv_path = CORRECTED_CSV if os.path.exists(CORRECTED_CSV) else BASE_CSV
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        total_blocks = 0
        reviewed_blocks = 0
        edit_ratios = []
        for label in BLOCK_LABELS:
            total_blocks += df[label].notna().sum()
            col = f'{label}_corrected'
            if col in df.columns:
                mask = df[col] == True
                reviewed_blocks += mask.sum()
                ratio_col = f'{label}_edit_ratio'
                if ratio_col in df.columns:
                    edit_ratios.extend(df.loc[mask, ratio_col].dropna().tolist())

        pct = 100 * reviewed_blocks / total_blocks if total_blocks else 0
        print(f'Progress: {reviewed_blocks} / {total_blocks} existing blocks '
              f'reviewed ({pct:.1f}%)')
        if edit_ratios:
            print(f'Mean edit ratio across reviewed blocks: {np.mean(edit_ratios):.3f}')

        bar_width = 50
        filled = int(bar_width * pct / 100)
        print('  [' + '#' * filled + '-' * (bar_width - filled) + f'] {pct:.1f}%')
    else:
        print('No results CSV found yet - run rsp_rate_by_block.py first.')

    print()
    print('=' * 78)


print('Starting live progress watcher. Press Ctrl+C to stop.\n')
try:
    while True:
        render_once()
        time.sleep(REFRESH_SECONDS)
except KeyboardInterrupt:
    print('\nStopped.')
