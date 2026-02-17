%% NEW_compare_units_by_correlation.m
% new script (re-structured by claude) using baseline normalized CCG values
% from NEW_make_ccgs.m to make summary figures for each correlation

% combines previous compare_units_by_correlation.m and make_summary_figs.m

%% compare_units_by_correlation.m
% Compare individual pairs of units - which units are coupled, do they share
% visual features, what is their cell identity, are receptive fields preserved
% across highly coupled pairs in wake, nrem, and rem sleep

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master")) 

sr = 30000;
nchan_probe = 384;
dtype = 'int16';

% kilosort directory
ksDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Elissa_Belluccini\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';

spikeTimes    = readNPY(fullfile(ksDir,'spike_times.npy'));
spikeTimes    = spikeTimes + 1;
spikeClusters = readNPY(fullfile(ksDir,'spike_clusters.npy'));
templates     = readNPY(fullfile(ksDir,'templates.npy'));
chanPos       = readNPY(fullfile(ksDir,'channel_positions.npy'));

cgFile = fullfile(ksDir, 'cluster_KSLabel.tsv');
cluster_groups = readtable(cgFile, 'FileType','text', 'Delimiter','\t');

keepGroups = {'good'};
toKeep     = ismember(cluster_groups.KSLabel, keepGroups);
keepClusters = cluster_groups.cluster_id(toKeep);

keepSpike     = ismember(spikeClusters, keepClusters);
spikeTimes    = spikeTimes(keepSpike);
spikeClusters = spikeClusters(keepSpike);

%% Get per unit depth
ypos = chanPos(:,2);

good_clusters = unique(spikeClusters);
nClusters     = length(good_clusters);

unitDepth = zeros(nClusters,1);
for i = 1:nClusters
    clu = good_clusters(i);
    tempIdx = mode(spikeClusters(spikeClusters==clu));
    template = squeeze(templates(tempIdx+1,:,:));
    [~,peakChan] = max(max(abs(template),[],2));
    unitDepth(i) = ypos(peakChan);
end

%% Load pre-computed CCGs and normalized coupling (from make_ccgs.m)
ccg_filePath = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\NEW_ccgs-all-units-40ms-Lag';
load(fullfile(ccg_filePath, 'ccgs_all'));
% This loads:
%   ccg_wake, ccg_nrem, ccg_rem, t          → short CCGs for visualization
%   ccg_wake_long, ccg_nrem_long, ccg_rem_long, t_long  → long CCGs
%   coupling_wake, coupling_nrem, coupling_rem           → normalized coupling
%   pairs, depthDiff, good_clusters, unitDepth           → metadata

states   = {'WAKE', 'NREM', 'REM'};
all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};
all_ccgs_long = {ccg_wake_long, ccg_nrem_long, ccg_rem_long};
% nBins    = size(ccg_wake, 1);
nPairs   = size(pairs, 1);

% Edge bins for visualization normalization (short CCG)
% nEdgeBins = floor(0.2 * nBins);

%% Sort pairs by normalized coupling strength
[sorted_coupling_wake, idx_wake] = sort(coupling_wake, 'descend');
[sorted_coupling_nrem, idx_nrem] = sort(coupling_nrem, 'descend');
[sorted_coupling_rem,  idx_rem]  = sort(coupling_rem,  'descend');

topN = 20;

%% Top coupled pairs during WAKE
disp('=== TOP COUPLED PAIRS DURING WAKE ===')
top_wake_table = table();
top_wake_table.Rank     = (1:topN)';
top_wake_table.UnitA    = pairs(idx_wake(1:topN), 1);
top_wake_table.UnitB    = pairs(idx_wake(1:topN), 2);
top_wake_table.ClusterA = good_clusters(pairs(idx_wake(1:topN), 1));
top_wake_table.ClusterB = good_clusters(pairs(idx_wake(1:topN), 2));
top_wake_table.Coupling = sorted_coupling_wake(1:topN);
top_wake_table.DepthA   = unitDepth(pairs(idx_wake(1:topN), 1));
top_wake_table.DepthB   = unitDepth(pairs(idx_wake(1:topN), 2));
top_wake_table.DepthDiff = abs(top_wake_table.DepthA - top_wake_table.DepthB);
disp(top_wake_table)

%% Top coupled pairs during NREM
disp('=== TOP COUPLED PAIRS DURING NREM ===')
top_nrem_table = table();
top_nrem_table.Rank     = (1:topN)';
top_nrem_table.UnitA    = pairs(idx_nrem(1:topN), 1);
top_nrem_table.UnitB    = pairs(idx_nrem(1:topN), 2);
top_nrem_table.ClusterA = good_clusters(pairs(idx_nrem(1:topN), 1));
top_nrem_table.ClusterB = good_clusters(pairs(idx_nrem(1:topN), 2));
top_nrem_table.Coupling = sorted_coupling_nrem(1:topN);
top_nrem_table.DepthA   = unitDepth(pairs(idx_nrem(1:topN), 1));
top_nrem_table.DepthB   = unitDepth(pairs(idx_nrem(1:topN), 2));
top_nrem_table.DepthDiff = abs(top_nrem_table.DepthA - top_nrem_table.DepthB);
disp(top_nrem_table)

%% Top coupled pairs during REM
disp('=== TOP COUPLED PAIRS DURING REM ===')
top_rem_table = table();
top_rem_table.Rank     = (1:topN)';
top_rem_table.UnitA    = pairs(idx_rem(1:topN), 1);
top_rem_table.UnitB    = pairs(idx_rem(1:topN), 2);
top_rem_table.ClusterA = good_clusters(pairs(idx_rem(1:topN), 1));
top_rem_table.ClusterB = good_clusters(pairs(idx_rem(1:topN), 2));
top_rem_table.Coupling = sorted_coupling_rem(1:topN);
top_rem_table.DepthA   = unitDepth(pairs(idx_rem(1:topN), 1));
top_rem_table.DepthB   = unitDepth(pairs(idx_rem(1:topN), 2));
top_rem_table.DepthDiff = abs(top_rem_table.DepthA - top_rem_table.DepthB);
disp(top_rem_table)

%% Find pairs most modulated between states
nrem_wake_modulation = coupling_nrem - coupling_wake;
rem_wake_modulation  = coupling_rem  - coupling_wake;

[~, idx_nrem_enhanced] = sort(nrem_wake_modulation, 'descend');

disp('=== PAIRS MOST ENHANCED DURING NREM (vs WAKE) ===')
nrem_enhanced_table = table();
nrem_enhanced_table.Rank          = (1:topN)';
nrem_enhanced_table.UnitA         = pairs(idx_nrem_enhanced(1:topN), 1);
nrem_enhanced_table.UnitB         = pairs(idx_nrem_enhanced(1:topN), 2);
nrem_enhanced_table.ClusterA      = good_clusters(pairs(idx_nrem_enhanced(1:topN), 1));
nrem_enhanced_table.ClusterB      = good_clusters(pairs(idx_nrem_enhanced(1:topN), 2));
nrem_enhanced_table.NREM_Coupling = coupling_nrem(idx_nrem_enhanced(1:topN));
nrem_enhanced_table.WAKE_Coupling = coupling_wake(idx_nrem_enhanced(1:topN));
nrem_enhanced_table.Modulation    = nrem_wake_modulation(idx_nrem_enhanced(1:topN));
nrem_enhanced_table.DepthDiff     = depthDiff(idx_nrem_enhanced(1:topN));
disp(nrem_enhanced_table)

%% Create pair summary PDFs
pdfDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\PAIR_SUMMARIES-TEST';

% Load additional data needed for the function
meanWaveformDir   = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters';
classificationDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SC-RGC-classification-with-kmeans-all-channels-BEST';

load(fullfile(meanWaveformDir,   'meanWav_units.mat'));
load(fullfile(classificationDir, 'RGC-SC-classification.mat'));

xpos = chanPos(:,1);
ypos = chanPos(:,2);

% Flip cell_type orientation if needed
cell_type  = cell_type';
cluster_id = cluster_id';

% Call the function
makePairSummaryPDFs( ...
    all_ccgs, ...
    all_ccgs_long, ...        % ADD - needed for normalization
    t, ...
    pairs, ...
    good_clusters, ...
    coupling_wake, ...
    coupling_nrem, ...
    coupling_rem, ...
    depthDiff, ...
    unitDepth, ...
    meanWav, ...
    cell_type, ...
    xpos, ...
    ypos, ...
    sr, ...
    baselineBins_long, ...    % ADD - needed for normalization
    binSize_short, ...        % ADD - needed for scaling
    binSize_long, ...         % ADD - needed for scaling
    'wake', ...
    20, ...
    pdfDir ...
);