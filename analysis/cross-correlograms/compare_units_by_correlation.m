%% Script to compare individual pairs of units
% - which units are coupled, do they share visual features, what is their
% cell identity, are receptive fields preserved across highly coupled pairs
% in wake, nrem, and rem sleep

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master")) 

sr = 30000;
nchan_probe = 384;
dtype = 'int16';

% kilosort directory
% ksDir = 'D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';
ksDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Elissa_Belluccini\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';

spikeTimes    = readNPY(fullfile(ksDir,'spike_times.npy'));
spikeTimes    = spikeTimes + 1;
spikeClusters = readNPY(fullfile(ksDir,'spike_clusters.npy')); % loads all units (not only the good ones)
templates     = readNPY(fullfile(ksDir,'templates.npy')); % can use template to compare? see if average waveforms from raw data are super different from  template
chanPos       = readNPY(fullfile(ksDir,'channel_positions.npy')); % [chan × 2]

% load cluster description (good/mua/noise)
% cgFile = fullfile(ksDir,'cluster_group.tsv'); % pre-manual curation
cgFile = fullfile(ksDir, 'cluster_KSLabel.tsv'); % post-manual curation

cluster_groups = readtable(cgFile, ...
    'FileType','text', ...
    'Delimiter','\t');

keepGroups = {'good'}; % only keep good labeled units

% if using post-manual cluster_KSLabel.tsv
toKeep = ismember(cluster_groups.KSLabel, keepGroups);
keepClusters = cluster_groups.cluster_id(toKeep); % contains only the good cluster ids (162 total)

% filter spikes by only kept clusters
keepSpike = ismember(spikeClusters, keepClusters);

spikeTimes    = spikeTimes(keepSpike);
spikeClusters = spikeClusters(keepSpike);

%% Get per unit depth
ypos = chanPos(:,2);

good_clusters = unique(spikeClusters);
nClusters = length(good_clusters);

unitDepth = zeros(nClusters,1);

for i = 1:nClusters
    clu = good_clusters(i);
    tempIdx = mode(spikeClusters(spikeClusters==clu));
    template = squeeze(templates(tempIdx+1,:,:)); % ks indexing
    
    [~,peakChan] = max(max(abs(template),[],2));
    unitDepth(i) = ypos(peakChan);
end

%% convert spikes to seconds
ts = double(spikeTimes) / sr;
gs = double(spikeClusters);

%% remap cluster ids into unit numbers
unitIDs = zeros(size(gs));
for i = 1:nClusters
    unitIDs(gs == good_clusters(i)) = i;
end

s = [ts unitIDs]; % to pass into CCG

%% Load Sleep States and Isolate by state
eegDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_eeg';
load(fullfile(eegDir, 'Mouse08_eeg.SleepState.states.mat'));

wakeInts = SleepState.ints.WAKEstate;
nremInts = SleepState.ints.NREMstate;
remInts  = SleepState.ints.REMstate;

s_wake = Restrict(s, wakeInts);
s_nrem = Restrict(s, nremInts);
s_rem  = Restrict(s, remInts);

% variables for ccg computation
binSize = 0.02; % s
duration = 8; % s

%% Load all CCGs (calculated in make_ccgs.m, saved in Mouse_08/ccgs-all-units

ccg_filePath = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\ccgs-all-units-40ms-Lag';
load(fullfile(ccg_filePath, 'ccgs_all'));

%% Build all the pairs - same as make_ccgs
pairs = nchoosek(1:nClusters,2);
nPairs = size(pairs,1);
depthDiff = zeros(nPairs,1);


%% extract coupling strength per pair of units based on depth - same as make_ccgs
zeroBin = ceil(length(t)/2);
win = zeroBin-2 : zeroBin+2;

coupling_wake = zeros(nPairs,1);
coupling_nrem = zeros(nPairs,1);
coupling_rem  = zeros(nPairs,1);

for p = 1:nPairs
    i = pairs(p,1);
    j = pairs(p,2);
    
    depthDiff(p) = abs(unitDepth(i) - unitDepth(j));
    
    coupling_wake(p) = sum(ccg_wake(win,i,j));
    coupling_nrem(p) = sum(ccg_nrem(win,i,j));
    coupling_rem(p)  = sum(ccg_rem(win,i,j));
end

% NEW
%% Find highest coupled pairs and get their cluster IDs 
%% THE CLUSTER ID IS THE RAW UNIT ID NUMBER USED IN MY VISUALIZATION SCRIPTS

% Sort pairs by coupling strength (choose which state to use)
[sorted_coupling_wake, idx_wake] = sort(coupling_wake, 'descend');
[sorted_coupling_nrem, idx_nrem] = sort(coupling_nrem, 'descend');
[sorted_coupling_rem, idx_rem] = sort(coupling_rem, 'descend');

% How many top pairs do you want to examine?
topN = 20;

%% Top coupled pairs during WAKE
disp('=== TOP COUPLED PAIRS DURING WAKE ===')
top_wake_table = table();
top_wake_table.Rank = (1:topN)';
top_wake_table.UnitA = pairs(idx_wake(1:topN), 1);
top_wake_table.UnitB = pairs(idx_wake(1:topN), 2);
top_wake_table.ClusterA = good_clusters(pairs(idx_wake(1:topN), 1));
top_wake_table.ClusterB = good_clusters(pairs(idx_wake(1:topN), 2));
top_wake_table.Coupling = sorted_coupling_wake(1:topN);
top_wake_table.DepthA = unitDepth(pairs(idx_wake(1:topN), 1));
top_wake_table.DepthB = unitDepth(pairs(idx_wake(1:topN), 2));
top_wake_table.DepthDiff = abs(top_wake_table.DepthA - top_wake_table.DepthB);

disp(top_wake_table)

%% Top coupled pairs during NREM
disp('=== TOP COUPLED PAIRS DURING NREM ===')
top_nrem_table = table();
top_nrem_table.Rank = (1:topN)';
top_nrem_table.UnitA = pairs(idx_nrem(1:topN), 1);
top_nrem_table.UnitB = pairs(idx_nrem(1:topN), 2);
top_nrem_table.ClusterA = good_clusters(pairs(idx_nrem(1:topN), 1));
top_nrem_table.ClusterB = good_clusters(pairs(idx_nrem(1:topN), 2));
top_nrem_table.Coupling = sorted_coupling_nrem(1:topN);
top_nrem_table.DepthA = unitDepth(pairs(idx_nrem(1:topN), 1));
top_nrem_table.DepthB = unitDepth(pairs(idx_nrem(1:topN), 2));
top_nrem_table.DepthDiff = abs(top_nrem_table.DepthA - top_nrem_table.DepthB);

disp(top_nrem_table)

%% Top coupled pairs during REM
disp('=== TOP COUPLED PAIRS DURING REM ===')
top_rem_table = table();
top_rem_table.Rank = (1:topN)';
top_rem_table.UnitA = pairs(idx_rem(1:topN), 1);
top_rem_table.UnitB = pairs(idx_rem(1:topN), 2);
top_rem_table.ClusterA = good_clusters(pairs(idx_rem(1:topN), 1));
top_rem_table.ClusterB = good_clusters(pairs(idx_rem(1:topN), 2));
top_rem_table.Coupling = sorted_coupling_rem(1:topN);
top_rem_table.DepthA = unitDepth(pairs(idx_rem(1:topN), 1));
top_rem_table.DepthB = unitDepth(pairs(idx_rem(1:topN), 2));
top_rem_table.DepthDiff = abs(top_rem_table.DepthA - top_rem_table.DepthB);

disp(top_rem_table)



















%% idk about this part
%% Plot CCGs for top N pairs in each state
figure;
for i = 1:min(topN, 9)  % plot up to 9 pairs
    subplot(3,3,i)
    
    % Get the pair indices
    unitA = pairs(idx_nrem(i), 1);
    unitB = pairs(idx_nrem(i), 2);
    clusterA = good_clusters(unitA);
    clusterB = good_clusters(unitB);
    
    % Plot CCGs for all states
    hold on;
    plot(t, squeeze(ccg_wake(:,unitA,unitB)), 'k', 'LineWidth', 1.5);
    plot(t, squeeze(ccg_nrem(:,unitA,unitB)), 'r', 'LineWidth', 1.5);
    plot(t, squeeze(ccg_rem(:,unitA,unitB)), 'b', 'LineWidth', 1.5);
    
    xlabel('Time lag (s)');
    ylabel('Spike count');
    title(sprintf('Pair #%d: Cluster %d & %d', i, clusterA, clusterB));
    legend({'WAKE', 'NREM', 'REM'}, 'Location', 'best');
    grid on;
    xlim([-0.5 0.5]);  % zoom into ±500 ms
end

sgtitle('Top Coupled Pairs (sorted by NREM coupling)')

%% Find pairs that are most modulated between states
% Calculate state modulation metrics
nrem_wake_modulation = coupling_nrem - coupling_wake;
rem_wake_modulation = coupling_rem - coupling_wake;
nrem_rem_modulation = coupling_nrem - coupling_rem;

% Find pairs with strongest NREM enhancement
[~, idx_nrem_enhanced] = sort(nrem_wake_modulation, 'descend');

disp('=== PAIRS MOST ENHANCED DURING NREM (vs WAKE) ===')
nrem_enhanced_table = table();
nrem_enhanced_table.Rank = (1:topN)';
nrem_enhanced_table.UnitA = pairs(idx_nrem_enhanced(1:topN), 1);
nrem_enhanced_table.UnitB = pairs(idx_nrem_enhanced(1:topN), 2);
nrem_enhanced_table.ClusterA = good_clusters(pairs(idx_nrem_enhanced(1:topN), 1));
nrem_enhanced_table.ClusterB = good_clusters(pairs(idx_nrem_enhanced(1:topN), 2));
nrem_enhanced_table.NREM_Coupling = coupling_nrem(idx_nrem_enhanced(1:topN));
nrem_enhanced_table.WAKE_Coupling = coupling_wake(idx_nrem_enhanced(1:topN));
nrem_enhanced_table.Modulation = nrem_wake_modulation(idx_nrem_enhanced(1:topN));
nrem_enhanced_table.DepthDiff = depthDiff(idx_nrem_enhanced(1:topN));

disp(nrem_enhanced_table)




