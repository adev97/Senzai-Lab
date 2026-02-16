%% Make cross-correlograms (CCGs) for different units to assess likelyhood of firing together

% **********************************
% calculates CCGs from CCG.m per pair per state (wake, nrem, rem).
% builds all pairs for the units
% creates heatmap with CCG score for all unique pairs
% includes short portion to plot CCGs for all states for one pair (manually
% assign two units to compare)


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
binSize = 0.001; % s (1ms bin for fine detail ccg indv plots) -- FOR HEATMAP - use 0.02s
duration = 0.040; % s (40ms lag for fine detail ccg indv plots) -- FOR HEATMAP - use 8s
% NOTE: this heatmap script does not work for small bin size + duration,
% does work for 4s instead of 8s though

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

% % Depth vs coupling
% figure
% scatter(depthDiff, coupling_nrem,'filled')
% xlabel('Depth difference (µm)')
% ylabel('NREM coupling strength')
% title('Coupling vs distance')
% 
% 
% % state modulation
% figure
% scatter(depthDiff, coupling_nrem - coupling_wake,'filled')
% xlabel('Depth difference (µm)')
% ylabel('NREM − WAKE coupling')
% title('Sleep-dependent coupling change')


%% 2. PLOT HEATMAP NEW BASELINE NORMALIZATION USE THIS ONE
%% Parameters
states = {'WAKE', 'NREM', 'REM'};
all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};
nStates = length(states);

nUnits = size(ccg_wake, 2); % number of units
nBins = size(ccg_wake, 1);

% Gaussian smoothing kernel (proper Gaussian)
sigma = 3; % in bins
xGauss = -3*sigma:3*sigma;
gaussKernel = exp(-xGauss.^2/(2*sigma^2));
gaussKernel = gaussKernel / sum(gaussKernel);

% Baseline window for normalization (e.g., 0.2-0.5 s away from zero lag)
baselineWindow = [0.2 0.5]; 

%% Build normalized pairwise CCGs (lower triangle only)
normCCG_states = cell(1, nStates);
pairCount = nUnits*(nUnits-1)/2;

for s = 1:nStates
    ccg = all_ccgs{s};
    
    normCCG = [];
    
    for i = 1:nUnits
        for j = 1:i-1  % lower triangle, each pair once
            pairCCG = squeeze(ccg(:,i,j))';
            
            % Smooth CCG with Gaussian kernel
            pairCCG_sm = conv(pairCCG, gaussKernel, 'same');
            
            % Find baseline bins
            baselineBins = find(abs(t) >= baselineWindow(1) & abs(t) <= baselineWindow(2));
            
            % Normalize by median of baseline
            pairCCGnorm = pairCCG_sm / median(pairCCG_sm(baselineBins));
            
            normCCG = [normCCG; pairCCGnorm];
        end
    end
    
    normCCG_states{s} = normCCG; % rows = pairs, columns = bins
end

%% Sort pairs by zero-lag coupling (from WAKE)
[~, zeroBin] = min(abs(t));  
zeroLagVals = normCCG_states{1}(:, zeroBin); 
[~, sortIdx] = sort(zeroLagVals, 'descend'); % highest coupling first

pairsToPlot = 400; % only top N pairs for readability % pairCount if want all
subsetIdx = sortIdx(1:pairsToPlot);

%% Plot heatmaps side-by-side
figure;

for s = 1:nStates
    subplot(1,nStates,s)
    
    dataSubset = normCCG_states{s}(subsetIdx,:);
    imagesc(t, 1:pairsToPlot, dataSubset);
    
    xlabel('Time lag (s)');
    ylabel('Neuron pairs (sorted by zero-lag WAKE)');
    title([states{s} ' normalized CCGs']);
    
    % Set color limits for consistent visualization
    clim([0.5 2.5]); 
    
    colormap('hot');
    colorbar;
    set(gca, 'YDir', 'reverse');
    xlim([t(1) t(end)]);
end


%% Plot CCG for one pair across all states with baseline normalization

% Specify the original cluster IDs you want to examine
clusterA = 63;  % original Kilosort cluster ID
clusterB = 415;  % original Kilosort cluster ID

% Convert to unit IDs
unitA = find(good_clusters == clusterA);
unitB = find(good_clusters == clusterB);

% Check if clusters exist in your good units
if isempty(unitA) || isempty(unitB)
    error('One or both cluster IDs not found in good_clusters');
end

disp(['Cluster ' num2str(clusterA) ' → Unit ' num2str(unitA)])
disp(['Cluster ' num2str(clusterB) ' → Unit ' num2str(unitB)])


figure; hold on;

colors = {'k','r','b'}; % color for each state
nEdgeBins = nBins;         % number of bins at each end to compute baseline

for s = 1:length(states)
    ccgState = all_ccgs{s};

    % extract CCG for this pair (spikes of B relative to A)
    ccgPair = squeeze(ccgState(:,unitA,unitB));

    % baseline normalization using edge bins
    baseline = median([ccgPair(1:nEdgeBins); ccgPair(end-nEdgeBins+1:end)]);
    ccgNorm = ccgPair / baseline;

    % plot
    plot(t, ccgNorm, 'Color', colors{s}, 'LineWidth', 2);
end

xlabel('Time lag (s)');
ylabel('Normalized CCG');
title(sprintf('Cluster %d → Cluster %d (Units %d → %d)', clusterA, clusterB, unitA, unitB));
legend(states, 'Location', 'best');
grid on;
xlim([t(1) t(end)]);


%% Save CCG calculations
ccgSavePath = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\ccgs-all-units-40ms-Lag";
if ~ccgSavePath
save(fullfile(ccgSavePath, 'ccgs_all.mat'), 'ccg_wake', 'ccg_nrem', 'ccg_rem', 't');
disp("CCGs (wake, nrem, rem) Saved!");









%% Example: Plot one unit pair across states with proper baseline and Gaussian smoothing

% % Choose a pair
% unit1 = 3; % first unit index
% unit2 = 7;  % second unit index
% 
% % Extract raw CCGs for this pair
% ccg_pair_wake = squeeze(ccg_wake(:, unit1, unit2));
% ccg_pair_nrem = squeeze(ccg_nrem(:, unit1, unit2));
% ccg_pair_rem  = squeeze(ccg_rem(:, unit1, unit2));
% 
% %% Define baseline bins (far from zero-lag, e.g., 200-500 ms)
% baselineBins = find(abs(t) > 0.2 & abs(t) < 0.5); % 200-500 ms away
% 
% %% Normalize by baseline median
% norm_pair_wake = ccg_pair_wake / median(ccg_pair_wake(baselineBins));
% norm_pair_nrem = ccg_pair_nrem / median(ccg_pair_nrem(baselineBins));
% norm_pair_rem  = ccg_pair_rem  / median(ccg_pair_rem(baselineBins));
% 
% %% Smooth with proper Gaussian
% sigma = 0.05; % 50 ms smoothing
% dt = t(2) - t(1); % bin width
% halfWidth = 4*sigma; 
% x = -halfWidth:dt:halfWidth;
% gaussKernel = exp(-x.^2/(2*sigma^2));
% gaussKernel = gaussKernel / sum(gaussKernel);
% 
% smoothCCG = @(c) conv(c, gaussKernel, 'same');
% 
% norm_pair_wake = smoothCCG(norm_pair_wake);
% norm_pair_nrem = smoothCCG(norm_pair_nrem);
% norm_pair_rem  = smoothCCG(norm_pair_rem);
% 
% %% Plot
% figure;
% 
% % Raw CCGs
% subplot(2,3,1)
% plot(t, ccg_pair_wake, 'b', 'LineWidth',1.5)
% xlabel('Time lag (s)'); ylabel('Raw count'); title('WAKE raw CCG')
% 
% subplot(2,3,2)
% plot(t, ccg_pair_nrem, 'r', 'LineWidth',1.5)
% xlabel('Time lag (s)'); ylabel('Raw count'); title('NREM raw CCG')
% 
% subplot(2,3,3)
% plot(t, ccg_pair_rem, 'g', 'LineWidth',1.5)
% xlabel('Time lag (s)'); ylabel('Raw count'); title('REM raw CCG')
% 
% % Normalized + smoothed CCGs
% subplot(2,3,4)
% plot(t, norm_pair_wake, 'b', 'LineWidth',1.5)
% xlabel('Time lag (s)'); ylabel('Normalized'); title('WAKE normalized')
% ylim([0 3]) % baseline=1, peaks visible
% 
% subplot(2,3,5)
% plot(t, norm_pair_nrem, 'r', 'LineWidth',1.5)
% xlabel('Time lag (s)'); ylabel('Normalized'); title('NREM normalized')
% ylim([0 3])
% 
% subplot(2,3,6)
% plot(t, norm_pair_rem, 'g', 'LineWidth',1.5)
% xlabel('Time lag (s)'); ylabel('Normalized'); title('REM normalized')
% ylim([0 3])


































































