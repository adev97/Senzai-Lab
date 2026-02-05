%% Make cross-correlograms (CCGs) for different units to assess likelyhood of firing together

% CCG.m from buzcode

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master")) 

% example function usage from yuta
% [ccg,t_ccg] = CCG(ts_res,gs_res,'binSize',0.02,'duration',8);

% goal: identify units that are firing together

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

keepGroups = {'good'}; % only keep good labeled units %%%%%%%%%%%%%%

% if using pre-manual cluster_group.tsv
% toKeep = ismember(cluster_groups.group, keepGroups);
% keepClusters = cluster_groups.cluster_id(toKeep); % contains only the good cluster ids (162 total)

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

% % convert s to [time unit]
% times = s.timestamps;      % cell array: {unit1, unit2, ...}
% uids  = s.units;        % unit IDs

% % Convert to FMAToolbox format
% allTimes = [];
% allIDs   = [];
% 
% for u = 1:length(times)
%     allTimes = [allTimes; times{u}(:)];
%     allIDs   = [allIDs;  u * ones(length(times{u}),1)];
% end
% 
% s_numeric = [allTimes allIDs];   % THIS is what Restrict wants


s_wake = Restrict(s, wakeInts);
s_nrem = Restrict(s, nremInts);
s_rem  = Restrict(s, remInts);

% variables for ccg computation
binSize = 0.02; % s
duration = 8; % s

%% Compute CCGs
[ccg_wake, t] = CCG(s_wake(:,1), s_wake(:,2), 'binSize',binSize,'duration',duration);
[ccg_nrem, ~] = CCG(s_nrem(:,1), s_nrem(:,2), 'binSize',binSize,'duration',duration);
[ccg_rem, ~] = CCG(s_rem(:,1), s_rem(:,2), 'binSize',binSize,'duration',duration);


%% Build all the pairs
pairs = nchoosek(1:nClusters,2);
nPairs = size(pairs,1);
depthDiff = zeros(nPairs,1);


%% extract coupling strength per pair of units based on depth
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


%% Visualize

% Depth vs coupling
figure
scatter(depthDiff, coupling_nrem,'filled')
xlabel('Depth difference (µm)')
ylabel('NREM coupling strength')
title('Coupling vs distance')


% state modulation
figure
scatter(depthDiff, coupling_nrem - coupling_wake,'filled')
xlabel('Depth difference (µm)')
ylabel('NREM − WAKE coupling')
title('Sleep-dependent coupling change')


%% Plot Heatmaps ordered by highest coupling - lowest coupling (from CCG peak)

%% Parameters
states = {'WAKE', 'NREM', 'REM'};
all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};
nStates = length(states);

nUnits = size(ccg_wake, 2); % number of units
nBins = size(ccg_wake, 1);

%% Build normalized pairwise CCGs (lower triangle only)
normCCG_states = cell(1, nStates);
pairCount = nUnits*(nUnits-1)/2;

for s = 1:nStates
    ccg = all_ccgs{s};
    
    normCCG = [];
    preID = [];
    postID = [];
    
    for i = 1:nUnits
        for j = 1:i-1  % lower triangle, each pair once
            pairCCG = squeeze(ccg(:,i,j))';
            pairCCGnorm = pairCCG / median(pairCCG); % normalize by median firing
            normCCG = [normCCG; pairCCGnorm];
            preID = [preID; i];
            postID = [postID; j];
        end
    end
    normCCG_states{s} = normCCG; % rows = pairs, columns = bins
end

%% Sort pairs by zero-lag coupling (from WAKE)
[~, zeroBin] = min(abs(t));  % bin closest to lag 0
zeroLagVals = normCCG_states{1}(:, zeroBin); % WAKE used for sorting
[~, sortIdx] = sort(zeroLagVals, 'descend'); % descending order (highest coupling to lowest coupling)

pairsToPlot = 400;
subsetIdx = sortIdx(1:pairsToPlot); % first N pairs after sorting

%% Plot heatmaps side-by-side
figure;

for s = 1:nStates
    subplot(1,nStates,s)

    dataSubset = normCCG_states{s}(subsetIdx, :);
    imagesc(t, 1:pairsToPlot, dataSubset);

    xlabel('Time lag (s)');
    ylabel('Neuron pairs (sorted by zero-lag WAKE)');
    title([states{s} ' normalized CCGs']);
    clim([0 2]);
    colormap('hot');
    colorbar;
    set(gca, 'YDir', 'reverse');
    xlim([t(1) t(end)]);
end


%% Plot CCG for one pair across all states
% choose your unit pair
unitA = 15;
unitB = 72;

figure; hold on;

colors = {'k','r','b'}; % color for each state

for s = 1:length(states)
    ccgState = all_ccgs{s};
    
    % extract lower triangle if needed
    ccgPair = squeeze(ccgState(:,unitA,unitB)); % spikes of B relative to A
    
    % normalize by median
    ccgNorm = ccgPair / median(ccgPair);
    
    plot(t, ccgNorm, 'Color', colors{s}, 'LineWidth', 2);
end

xlabel('Time lag (s)');
ylabel('Normalized CCG');
title(sprintf('Unit %d → Unit %d across states', unitA, unitB));
legend(states,'Location','best');
grid on;
xlim([t(1) t(end)]);




































































% fill x position and y position per channel
xpos = chanPos(:,1);
ypos = chanPos(:,2);

good_clusters = unique(spikeClusters);
good_clusters(good_clusters==0) = [];   % remove noise if present
nClusters = numel(good_clusters);

%% spikeTimes - contains all spike timings (when)
%% spikeClusters - contains which unit gave the spike (who) and preserves channel number structure

%% Now building cross-correlograms, searching for units that fire together
% during this data collection

% convert time to seconds instead of samples
ts = double(spikeTimes) / sr;   % convert to seconds
gs = double(spikeClusters);     % unit IDs

%% remap cluster ids into unit numbers
unitIDs = zeros(size(gs));
for i = 1:nClusters
    unitIDs(gs == good_clusters(i)) = i;
end

s = [ts unitIDs]; % to pass into CCG

% % Map original cluster IDs to consecutive indices
% [~, gs_idx] = ismember(spikeClusters, good_clusters);
% 
% % % Use gs_idx for CCG
% % [ccg, t] = CCG(ts, gs_idx, 'binSize', binSize, 'duration', duration);

%% Open Sleep States files

eegDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_eeg';
load(fullfile(eegDir, 'Mouse08_eeg.SleepState.states.mat'));

% set bin size and duration to find units firing together
binSize  = 0.01;   % 10 ms resolution (default in CCG.m)
duration = 8;    % ±8 ms window (same as yuta's example script)

%% Isolate different sleep states
s_wake = Restrict(s, SleepState.ints.WAKEstate);
s_nrem = Restrict(s, SleepState.ints.NREMstate);
s_rem  = Restrict(s, SleepState.ints.REMstate);

%% Compute CCGs for all states
states = {'WAKEstate','REMstate','NREMstate'};

for st = 1:length(states)

    intervals = SleepState.ints.(states{st});

    % keep only spikes in that state
    s_state = Restrict(s, intervals);

    % compute cross-correlogram
    [ccg, t] = CCG(s_state(:,1), s_state(:,2), ...
                   'binSize', binSize, ...
                   'duration', duration);

    CCG_all.(states{st}).ccg = ccg;
    CCG_all.(states{st}).t   = t;
end

%% Make heatmap
for st = 1:length(states)

    ccg = CCG_all.(states{st}).ccg;
    t   = CCG_all.(states{st}).t;

    CCG_pooled = [];

    for k = 1:nClusters
        CCGmtrx = squeeze(ccg(:,k,1:(k-1)))';   % avoid duplicates
        CCG_pooled = [CCG_pooled; CCGmtrx];
    end

    % Normalize (reveals real coupling)
    CCG_pooled = CCG_pooled ./ median(CCG_pooled,2);

    figure;
    imagesc(t*1000, 1:size(CCG_pooled,1), CCG_pooled);
    xlabel('Lag (ms)');
    ylabel('Neuron pairs');
    title(states{st});
    colorbar;
end



