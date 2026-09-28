% rsp_breathmetrics_by_condition.m
%
% Respiratory rate per participant, per block, per condition, computed
% with BreathMetrics (Noto et al., 2018) instead of NeuroKit2 - run on
% the WHOLE continuous recording (not pre-cut segments), then breaths are
% assigned to blocks afterward using the trigger markers. This matters:
% BreathMetrics' baseline correction and multi-window peak detection
% benefit from seeing the full continuous signal, not an artificially
% truncated piece - cutting first can distort the baseline/edges.
%
% Data source: .mat files exported by export_rsp_for_matlab.py, containing:
%   rsp            - raw respiration signal (Sensor-H:RSP)
%   sfreq          - sampling rate (Hz)
%   marker_onset   - onset time (s) of every trigger marker
%   marker_desc    - text of every trigger marker, e.g. "1,1" (block,condition)
%   row, name, started_with - identity/behavioral-list info
%
% Marker code convention (block-first, confirmed earlier in this project):
%   first digit = block (1 or 2), second digit = condition (1=certain, 2=uncertain)
%   each block is marked by the SAME code twice: once at start, once at end.
%
% Sensor type: chestbelt (chest expansion), so BreathMetrics dataType =
% 'humanBB' (breathing belt), NOT 'humanAirflow' (that's for nasal
% cannula/airflow sensors - a different signal morphology).

clear; clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));

BREATHMETRICS_ROOT = fullfile(SCRIPT_DIR, 'breathmetrics-master');
MAT_DIR = fullfile(SCRIPT_DIR, 'data', 'rsp_mat_export');
OUT_DIR = fullfile(SCRIPT_DIR, 'output', 'rsp_breathmetrics_output');

% Machine-specific overrides, gitignored - copy local_config.m.example to
% local_config.m to point these at your own data/toolbox locations
% instead of editing this file.
LOCAL_CONFIG = fullfile(SCRIPT_DIR, 'local_config.m');
if exist(LOCAL_CONFIG, 'file')
    run(LOCAL_CONFIG);
end

addpath(genpath(BREATHMETRICS_ROOT));
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end
if ~exist(fullfile(OUT_DIR, 'diagnostic_plots'), 'dir')
    mkdir(fullfile(OUT_DIR, 'diagnostic_plots'));
end

DATA_TYPE = 'humanBB';  % chestbelt signal, not nasal airflow

matFiles = dir(fullfile(MAT_DIR, '*.mat'));
fprintf('Found %d subjects to process\n\n', numel(matFiles));

results = struct('row', {}, 'name', {}, 'started_with', {}, ...
    'rate_certain', {}, 'rate_uncertain', {}, ...
    'rate_certain_1', {}, 'rate_certain_2', {}, ...
    'rate_uncertain_1', {}, 'rate_uncertain_2', {}, ...
    'nBreaths_certain', {}, 'nBreaths_uncertain', {}, ...
    'flatline_pct', {}, 'clip_pct', {}, 'pause_pct_overall', {});

for fi = 1:numel(matFiles)
    matPath = fullfile(matFiles(fi).folder, matFiles(fi).name);
    D = load(matPath);

    rsp = double(D.rsp(:)');
    sfreq = double(D.sfreq);
    markerOnset = double(D.marker_onset(:));
    markerDesc = cellstr(D.marker_desc);
    rowNum = double(D.row);
    name = char(D.name);
    startedWith = double(D.started_with);

    fprintf('--- %02d %s ---\n', rowNum, name);

    % -------------------------------------------------------------
    % Basic signal-quality checks (same logic as the earlier Python
    % quality report: flatline = stuck sensor, clipping = saturation)
    % -------------------------------------------------------------
    diffs = diff(rsp);
    flatlinePct = sum(diffs == 0) / numel(diffs) * 100;
    rng = max(rsp) - min(rsp);
    tol = rng * 0.001;
    clipPct = (sum(rsp >= max(rsp) - tol) + sum(rsp <= min(rsp) + tol)) / numel(rsp) * 100;

    % -------------------------------------------------------------
    % Run BreathMetrics on the FULL continuous recording
    % -------------------------------------------------------------
    try
        bmObj = breathmetrics(rsp, sfreq, DATA_TYPE);
        bmObj.estimateAllFeatures();
    catch err
        fprintf('  BreathMetrics FAILED: %s\n', err.message);
        continue
    end

    % inhale onset times (seconds) - one per detected breath, used to
    % assign each breath to whichever block/condition window it falls in
    inhaleOnsetTimes = bmObj.time(bmObj.inhaleOnsets);

    secondaryFeatures = bmObj.secondaryFeatures;
    if isKey(secondaryFeatures, 'Percent Of Breaths With Inhale Pause')
        pausePct = secondaryFeatures('Percent Of Breaths With Inhale Pause');
    else
        pausePct = NaN;
    end

    % -------------------------------------------------------------
    % Parse markers into block segments: {block,cond} -> [start_s end_s]
    % -------------------------------------------------------------
    [~, order] = sort(markerOnset);
    markerOnset = markerOnset(order);
    markerDesc = markerDesc(order);

    pendingKeys = {};
    pendingTimes = [];
    segCertain = [];   % Nx2 [start end]
    segUncertain = [];

    for k = 1:numel(markerDesc)
        tok = regexp(markerDesc{k}, '^([12]),([12])$', 'tokens');
        if isempty(tok)
            continue
        end
        block = tok{1}{1};
        cond = tok{1}{2};
        key = [block '_' cond];
        idxExisting = find(strcmp(pendingKeys, key), 1);
        if isempty(idxExisting)
            pendingKeys{end+1} = key; %#ok<AGROW>
            pendingTimes(end+1) = markerOnset(k); %#ok<AGROW>
        else
            startT = pendingTimes(idxExisting);
            endT = markerOnset(k);
            pendingKeys(idxExisting) = [];
            pendingTimes(idxExisting) = [];
            if strcmp(cond, '1')
                segCertain = [segCertain; startT endT]; %#ok<AGROW>
            else
                segUncertain = [segUncertain; startT endT]; %#ok<AGROW>
            end
        end
    end

    % -------------------------------------------------------------
    % Breathing rate per condition: count breaths (inhale onsets) whose
    % time falls inside that condition's block windows, divide by total
    % time covered by those blocks -> breaths per minute.
    % -------------------------------------------------------------
    [rateCertain, nCertain] = rateFromSegments(inhaleOnsetTimes, segCertain);
    [rateUncertain, nUncertain] = rateFromSegments(inhaleOnsetTimes, segUncertain);

    % Per-block rates too (not just the 2-block average): segCertain/segUncertain
    % are built in chronological order of completion, and since block 1's two
    % segments (certain+uncertain) both finish before block 2 starts, row 1 is
    % always block 1's segment and row 2 is always block 2's segment.
    rateCertain1 = NaN; rateCertain2 = NaN;
    rateUncertain1 = NaN; rateUncertain2 = NaN;
    if size(segCertain, 1) >= 1
        rateCertain1 = rateFromSegments(inhaleOnsetTimes, segCertain(1,:));
    end
    if size(segCertain, 1) >= 2
        rateCertain2 = rateFromSegments(inhaleOnsetTimes, segCertain(2,:));
    end
    if size(segUncertain, 1) >= 1
        rateUncertain1 = rateFromSegments(inhaleOnsetTimes, segUncertain(1,:));
    end
    if size(segUncertain, 1) >= 2
        rateUncertain2 = rateFromSegments(inhaleOnsetTimes, segUncertain(2,:));
    end

    fprintf('  certain: rate=%.1f bpm (n=%d breaths)   uncertain: rate=%.1f bpm (n=%d breaths)\n', ...
        rateCertain, nCertain, rateUncertain, nUncertain);
    fprintf('  block1: certain=%.1f uncertain=%.1f | block2: certain=%.1f uncertain=%.1f\n', ...
        rateCertain1, rateUncertain1, rateCertain2, rateUncertain2);
    fprintf('  flatline=%.2f%%  clip=%.2f%%  pause%%=%.1f\n', flatlinePct, clipPct, pausePct);

    results(end+1) = struct('row', rowNum, 'name', name, 'started_with', startedWith, ...
        'rate_certain', rateCertain, 'rate_uncertain', rateUncertain, ...
        'rate_certain_1', rateCertain1, 'rate_certain_2', rateCertain2, ...
        'rate_uncertain_1', rateUncertain1, 'rate_uncertain_2', rateUncertain2, ...
        'nBreaths_certain', nCertain, 'nBreaths_uncertain', nUncertain, ...
        'flatline_pct', flatlinePct, 'clip_pct', clipPct, ...
        'pause_pct_overall', pausePct); %#ok<AGROW>

    % -------------------------------------------------------------
    % Diagnostic plot: full signal, detected inhale peaks, block
    % boundaries shaded - so each subject can be visually spot-checked
    % without re-running BreathMetrics.
    % -------------------------------------------------------------
    fig = figure('Visible', 'off', 'Position', [100 100 1400 500]);
    hold on;
    plot(bmObj.time, bmObj.baselineCorrectedRespiration, 'k-', 'LineWidth', 0.5);
    scatter(bmObj.time(bmObj.inhalePeaks), bmObj.peakInspiratoryFlows, 20, 'r', 'filled');
    for r = 1:size(segCertain, 1)
        xline(segCertain(r,1), 'b--'); xline(segCertain(r,2), 'b--');
    end
    for r = 1:size(segUncertain, 1)
        xline(segUncertain(r,1), 'r--'); xline(segUncertain(r,2), 'r--');
    end
    title(sprintf('%02d %s (blue=certain blocks, red=uncertain blocks)', rowNum, name), 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('Baseline-corrected RSP');
    saveas(fig, fullfile(OUT_DIR, 'diagnostic_plots', sprintf('%02d_%s.png', rowNum, name)));
    close(fig);
end

% ---------------------------------------------------------------------
% Save summary table
% ---------------------------------------------------------------------
T = struct2table(results);
writetable(T, fullfile(OUT_DIR, 'breathmetrics_rate_by_condition.csv'));
fprintf('\nSaved summary: %s\n', fullfile(OUT_DIR, 'breathmetrics_rate_by_condition.csv'));
fprintf('Saved %d diagnostic plots to: %s\n', numel(matFiles), fullfile(OUT_DIR, 'diagnostic_plots'));

% ---------------------------------------------------------------------
% Flag any subject/block whose rate falls outside a plausible 10-18
% breaths/min range - checked per INDIVIDUAL block, not the 2-block
% condition average, so a bad block doesn't get hidden by averaging.
% ---------------------------------------------------------------------
RATE_LOW = 10;
RATE_HIGH = 18;
blockCols = {'rate_certain_1', 'rate_uncertain_1', 'rate_certain_2', 'rate_uncertain_2'};
blockNames = {'1st Certain', '1st Uncertain', '2nd Certain', '2nd Uncertain'};

fprintf('\n=== Blocks outside %g-%g breaths/min ===\n', RATE_LOW, RATE_HIGH);
nOutOfRange = 0;
for i = 1:height(T)
    for b = 1:4
        val = T.(blockCols{b})(i);
        if ~isnan(val) && (val < RATE_LOW || val > RATE_HIGH)
            fprintf('%02d-%-16s %-14s rate=%.1f bpm\n', T.row(i), T.name{i}, blockNames{b}, val);
            nOutOfRange = nOutOfRange + 1;
        end
    end
end
fprintf('Total out-of-range blocks: %d (out of %d subjects x 4 blocks = %d)\n', ...
    nOutOfRange, height(T), height(T) * 4);

% ---------------------------------------------------------------------
% Violin + dots plot, same style as the earlier PhysioNet HRV scripts,
% with each dot labeled by its subject ID (row-name).
% ---------------------------------------------------------------------
figure('Position', [100 100 800 800]);
hold on;
cVals = T.rate_certain;
uVals = T.rate_uncertain;
plotViolin(1, cVals(~isnan(cVals)), [0.75 0.75 0.75]);
plotViolin(2, uVals(~isnan(uVals)), [0.75 0.75 0.75]);
for i = 1:height(T)
    if isnan(cVals(i)) || isnan(uVals(i))
        continue
    end
    if T.started_with(i) == 1
        col = [0 0.3 1];
    elseif T.started_with(i) == 2
        col = [1 0.1 0.1];
    else
        col = [0.6 0.6 0.6];
    end
    label = sprintf('%02d-%s', T.row(i), T.name{i});
    plot([1 2], [cVals(i) uVals(i)], '-', 'Color', [col 0.25], 'LineWidth', 0.8);
    scatter(1, cVals(i), 40, col, 'filled');
    scatter(2, uVals(i), 40, col, 'filled');
    text(1 - 0.06, cVals(i), label, 'FontSize', 6, 'HorizontalAlignment', 'right');
    text(2 + 0.06, uVals(i), label, 'FontSize', 6, 'HorizontalAlignment', 'left');
end
xlim([0.3 2.7]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Certain', 'Uncertain'});
ylabel('Breathing rate (breaths/min), BreathMetrics');
title('Respiratory rate by condition - BreathMetrics (humanBB)');
grid on;
hold off;

% ---------------------------------------------------------------------
% Block-by-block plot: 1st certain, 1st uncertain, 2nd certain, 2nd
% uncertain, in that chronological/labeled order - shows the within-
% subject progression across all 4 individual blocks, not just the
% 2-block-averaged condition comparison above.
% ---------------------------------------------------------------------
figure('Position', [100 100 1000 800]);
hold on;
blockVals = {T.rate_certain_1, T.rate_uncertain_1, T.rate_certain_2, T.rate_uncertain_2};
blockLabels = {'1st Certain', '1st Uncertain', '2nd Certain', '2nd Uncertain'};
for b = 1:4
    v = blockVals{b};
    plotViolin(b, v(~isnan(v)), [0.75 0.75 0.75]);
end
for i = 1:height(T)
    vals = [T.rate_certain_1(i), T.rate_uncertain_1(i), T.rate_certain_2(i), T.rate_uncertain_2(i)];
    if any(isnan(vals))
        continue
    end
    if T.started_with(i) == 1
        col = [0 0.3 1];
    elseif T.started_with(i) == 2
        col = [1 0.1 0.1];
    else
        col = [0.6 0.6 0.6];
    end
    label = sprintf('%02d-%s', T.row(i), T.name{i});
    plot(1:4, vals, '-', 'Color', [col 0.2], 'LineWidth', 0.8);
    scatter(1:4, vals, 35, col, 'filled');
    text(1 - 0.08, vals(1), label, 'FontSize', 6, 'HorizontalAlignment', 'right');
    text(4 + 0.08, vals(4), label, 'FontSize', 6, 'HorizontalAlignment', 'left');
end
xlim([0.4 4.6]);
set(gca, 'XTick', 1:4, 'XTickLabel', blockLabels);
ylabel('Breathing rate (breaths/min), BreathMetrics');
title('Respiratory rate by block (1st certain -> 1st uncertain -> 2nd certain -> 2nd uncertain)');
grid on;
hold off;

% rateFromSegments() and plotViolin() are now separate files in this same
% folder (rateFromSegments.m, plotViolin.m) rather
% than local functions at the end of this script - local functions in a
% script only get registered if the ENTIRE file is run as one unit, and
% separate files avoid that failure mode entirely.
