"""
Shared marker-parsing logic, imported by rsp_rate_by_block.py,
manual_correct_single.py, and manual_correct_peaks.py, so all three
scripts always agree on block boundaries.

Two layers:
1. GENERAL parser: accepts both comma ("1,2") and period ("1.2") coded
   markers, debounces double-tapped markers (repeat within 10s of the
   same code = a bounce, not a real end), and falls back to a
   single-occurrence convention (block ends at the next marker in the
   file) for codes that never get a matching second occurrence -
   gated by a 60s minimum plausible duration so a bad fallback pairing
   is discarded rather than kept as a garbage micro-block.
2. PER-SUBJECT OVERRIDES: a small number of subjects have marker
   mistakes that no general rule can safely infer (a marker recorded
   with the wrong code, or a marker missing entirely). These were
   diagnosed by hand by inspecting each subject's raw marker sequence
   and are applied on top of the general parser's output, replacing or
   filling in the specific blocks affected. Every override below is
   documented with exactly what was wrong and the reasoning for the fix.
   Overrides are keyed by zero-padded subject ID (matched against the
   leading "<ID>_" of the EDF filename), not by participant name, so
   this file carries no identifying information.
"""
import re

CODE_RE = re.compile(r'^([12])[,.]([12])$')
BOUNCE_TOLERANCE_S = 10
MIN_FALLBACK_DURATION_S = 60
SUBJECT_ID_RE = re.compile(r'^(\d+)_')


def _general_parse(raw):
    onsets = raw.annotations.onset
    descs = raw.annotations.description
    order = onsets.argsort()
    onsets = onsets[order]
    descs = [descs[i] for i in order]

    pending = {}
    segments = {}
    used_idx = set()
    for i, (onset, desc) in enumerate(zip(onsets, descs)):
        m = CODE_RE.match(desc)
        if not m:
            continue
        key = f'{m.group(1)}_{m.group(2)}'
        if key not in pending:
            pending[key] = (onset, i)
        else:
            start_t, start_i = pending[key]
            if onset - start_t <= BOUNCE_TOLERANCE_S:
                pending[key] = (onset, i)  # debounce
                continue
            segments[key] = (start_t, onset)
            used_idx.add(start_i)
            used_idx.add(i)
            del pending[key]

    still_pending = {k: v for k, v in pending.items() if k not in segments}
    for i, (onset, desc) in enumerate(zip(onsets, descs)):
        if i in used_idx:
            continue
        m = CODE_RE.match(desc)
        if not m:
            continue
        key = f'{m.group(1)}_{m.group(2)}'
        if key in segments:
            continue
        if key in still_pending and still_pending[key][1] != i:
            continue
        if i + 1 < len(onsets):
            end = onsets[i + 1]
            if end - onset >= MIN_FALLBACK_DURATION_S:
                segments[key] = (onset, end)
        still_pending.pop(key, None)

    return segments


# ---------------------------------------------------------------------------
# Per-subject overrides. Each entry replaces/fills specific block keys
# for that exact filename. Diagnosed by hand from the raw marker text -
# see the comment above each entry for what was actually wrong.
# ---------------------------------------------------------------------------
MANUAL_OVERRIDES = {
    # Markers use bare '1'/'2'/'1y' with no comma/period at all - the
    # general parser can't recognize these as codes. True structure
    # (confirmed by hand): '1y' (debounced start) -> '1' = 1st Certain;
    # '2' -> '2' = 1st Uncertain; '1' -> '1' = 2nd Certain;
    # '2' -> '1' = 2nd Uncertain (last block, uses whatever markers remain).
    '01': {
        '1_1': (214.33, 442.58),
        '1_2': (450.63, 661.84),
        '2_1': (669.77, 872.23),
        '2_2': (880.79, 1077.23),
    },
    # Confirmed true order: 1,2...1,2 (1st Uncertain) -> 1,1...1,1 (1st
    # Certain) -> 2,1...2,1 (2nd Certain, ALREADY correctly resolved by
    # the general parser from the two untouched real '2,1' markers - not
    # overridden here, so the existing saved correction for this block
    # stays intact) -> 2,2...2,2 (2nd Uncertain). The marker at t=700.80
    # was recorded as '1,2' but should have been '2,2', pairing with the
    # existing '2,2'@870.70 to recover 2nd Uncertain, which was
    # previously unresolved (only one usable raw marker).
    '05': {
        '2_2': (700.80, 870.70),
    },
    # 2nd Certain's start marker was never logged - only its end
    # ('2,1'@896.32) exists. Estimated start using the ~10s gap this
    # same file shows elsewhere between a block's end and the next
    # block's start (449.55 -> 459.27 = 9.72s), applied to the end of
    # 2nd Uncertain (754.03 + ~10s).
    '15': {
        '2_1': (764.03, 896.32),
    },
    # The marker at t=245.91 was recorded as '2,1' but should have been
    # '1,2' - true block order is 1,1...1,1 then 1,2...1,2 (both block-1
    # conditions complete before block 2 starts). Relabeling this one
    # marker lets the general parser's own pairing/fallback correctly
    # recover both 1st Uncertain and 2nd Certain.
    '26': {
        '1_2': (245.91, 393.76),
        '2_1': (403.44, 583.03),
    },
    # The marker at t=891.06 was recorded as '2,1' but should have been
    # '1,1', pairing with the existing '1,1'@673.48 to close 1st Certain.
    # IMPORTANT: because the general parser only sees the raw ORIGINAL
    # text, it independently (and wrongly) paired that same 891.06
    # marker with the next '2,1'@1105.37 to form a stale "2nd Certain"
    # (891.06-1105.37) - that has to be explicitly overridden here too,
    # not just '1_1', or the stale value survives untouched. True 2nd
    # Certain is 1105.37-1278.91 (the marker after 'Artifact!', which is
    # not a code marker and is correctly ignored either way). 2nd
    # Uncertain (898.50-1101.94) was already correct and is unaffected.
    '39': {
        '1_1': (673.48, 891.06),
        '2_1': (1105.37, 1278.91),
    },
    # Both 1st Uncertain's start and 2nd Uncertain's start were never
    # logged - only their end markers exist ('1,2'@696.98 and
    # '2,2'@964.76). Estimated each start using the same ~10s gap
    # convention as subject 15, applied right after the preceding
    # block's end (1st Certain ends 504.91; 2nd Certain ends 817.66).
    '42': {
        '1_2': (514.91, 696.98),
        '2_2': (827.66, 964.76),
    },
    # 2nd Certain's start was never logged - only its end
    # ('2,1'@939.59) exists. Estimated start using the ~10s gap
    # convention, applied after 1st Uncertain ends (745.87).
    '43': {
        '2_1': (755.87, 939.59),
    },
    # The marker at t=277.22 is a stray double-tap of '1,1' just outside
    # the general parser's normal 10s debounce window (10.53s later),
    # not a real block end. True 1st Certain pairs the first '1,1'
    # occurrence (266.69) with the third one (545.42); the erroneous
    # middle occurrence is ignored entirely. 1st Uncertain, 2nd Certain,
    # and 2nd Uncertain were already correctly resolved and are
    # unaffected.
    '45': {
        '1_1': (266.69, 545.42),
    },
    # 1st Uncertain's start was never logged - only its end
    # ('1,2'@734.77) exists. Estimated start using the ~10s gap
    # convention, applied after 1st Certain ends (512.51).
    '47': {
        '1_2': (522.51, 734.77),
    },
    # 1st Certain's start was never logged - only its end ('1,1'@428.74)
    # exists. Estimated start using the ~10s gap convention, applied
    # after 1st Uncertain ends (297.73).
    '48': {
        '1_1': (307.73, 428.74),
    },
    # The marker at t=491.33 ("the first 2,1") was recorded as '2,1' but
    # should have been '1,2', pairing with the existing '1,2'@300.02 to
    # close 1st Uncertain. IMPORTANT: the general parser's own default
    # (using the raw original text) paired that same 491.33 marker with
    # the very next '2,1'@502.15, producing a stale, bogus ~11s "2nd
    # Certain" - that has to be explicitly overridden too. True 2nd
    # Certain pairs '2,1'@502.15 with '2,1'@704.25. 1st Certain
    # (60.92-290.11) and 2nd Uncertain (714.38-882.97) were already
    # correct and are unaffected.
    '51': {
        '1_2': (300.02, 491.33),
        '2_1': (502.15, 704.25),
    },
    # The marker at t=1347.84 was recorded as '2,2' but should have been
    # '2,1', pairing with the existing '2,1'@1628.06 to close 2nd
    # Certain (previously unresolved - only one usable raw '2,1'
    # marker). That leaves the two remaining real '2,2' markers
    # (1638.77, 1888.50) to pair naturally as 2nd Uncertain, with no
    # estimation needed. 1st Certain (529.88-1011.90) and 1st Uncertain
    # (1020.12-1337.75) were already correct and are unaffected. The
    # several 'Artifact!'/'artifact'/'art' markers scattered through
    # this file are not code markers and are correctly ignored either way.
    '52': {
        '2_1': (1347.84, 1628.06),
        '2_2': (1638.77, 1888.50),
    },
}


def get_block_segments(raw, filename=None):
    """
    Returns {'1_1': (start_s, end_s), ...} for whichever of the 4 blocks
    could be resolved. filename (just the basename, e.g. "05_subjectname.edf")
    is used to apply hand-diagnosed per-subject overrides on top of the
    general parser's output - pass it whenever available. Overrides are
    matched by the leading numeric subject ID, not the full filename.
    """
    segments = _general_parse(raw)
    if filename:
        m = SUBJECT_ID_RE.match(filename)
        if m and m.group(1) in MANUAL_OVERRIDES:
            segments.update(MANUAL_OVERRIDES[m.group(1)])
    return segments
