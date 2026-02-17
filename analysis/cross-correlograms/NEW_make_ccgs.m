%% NEW_make_ccgs.m
% new script (restructured by claude) using baseline normalized CCG values
% based on make_ccgs.m

% calulcates a long window from which to get baseline coincidence firing

%% make_ccgs.m
% Make cross-correlograms (CCGs) for different units to assess likelihood of firing together

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
chanPos       = readNPY(fullfile(ksDir,'channel_positions.npy')); % [chan x 2]

cgFile = fullfile(ksDir, 'cluster_KSLabel.tsv');
cluster_groups = readtable(cgFile, 'FileType','text', 'Delimiter','\t');

keepGroups = {'good'};
toKeep = ismember(cluster_groups.KSLabel, keepGroups);
keepClusters = cluster_groups.cluster_id(toKeep);

keepSpike = ismember(spikeClusters, keepClusters);
spikeTimes    = spikeTimes(keepSpike);
spikeClusters = spikeClusters(keepSpike);

%% Get per unit depth
xpos = chanPos(:,1);
ypos = chanPos(:,2);

good_clusters = unique(spikeClusters);
nClusters = length(good_clusters);

unitDepth = zeros(nClusters,1);

for i = 1:nClusters
    clu = good_clusters(i);
    tempIdx = mode(spikeClusters(spikeClusters==clu));
    template = squeeze(templates(tempIdx+1,:,:));
    [~,peakChan] = max(max(abs(template),[],2));
    unitDepth(i) = ypos(peakChan);
end

%% Convert spikes to seconds and remap cluster IDs to unit numbers
ts = double(spikeTimes) / sr;
gs = double(spikeClusters);

unitIDs = zeros(size(gs));
for i = 1:nClusters
    unitIDs(gs == good_clusters(i)) = i;
end

s = [ts unitIDs];

%% Load Sleep States
eegDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_eeg';
load(fullfile(eegDir, 'Mouse08_eeg.SleepState.states.mat'));

wakeInts = SleepState.ints.WAKEstate;
nremInts = SleepState.ints.NREMstate;
remInts  = SleepState.ints.REMstate;

s_wake = Restrict(s, wakeInts);
s_nrem = Restrict(s, nremInts);
s_rem  = Restrict(s, remInts);

%% ===== SHORT CCG (for visualization) =====
% 1ms bins, 40ms window - used for plotting individual pairs
binSize_short  = 0.001; % 1ms bins
duration_short = 0.040; % 40ms total window

fprintf('Computing short CCGs for visualization...\n');
[ccg_wake, t] = CCG(s_wake(:,1), s_wake(:,2), 'binSize', binSize_short, 'duration', duration_short);
[ccg_nrem, ~] = CCG(s_nrem(:,1), s_nrem(:,2), 'binSize', binSize_short, 'duration', duration_short);
[ccg_rem,  ~] = CCG(s_rem(:,1),  s_rem(:,2),  'binSize', binSize_short, 'duration', duration_short);
fprintf('Short CCGs done.\n');

%% ===== LONG CCG (for normalized coupling calculation) =====
% 1ms bins, 500ms window - gives far-away bins for valid baseline
binSize_long  = 0.020; % 20ms bins
duration_long = 8; % 8s window - far enough for baseline

fprintf('Computing long CCGs for normalized coupling...\n');
[ccg_wake_long, t_long] = CCG(s_wake(:,1), s_wake(:,2), 'binSize', binSize_long, 'duration', duration_long);
[ccg_nrem_long, ~]      = CCG(s_nrem(:,1), s_nrem(:,2), 'binSize', binSize_long, 'duration', duration_long);
[ccg_rem_long,  ~]      = CCG(s_rem(:,1),  s_rem(:,2),  'binSize', binSize_long, 'duration', duration_long);
fprintf('Long CCGs done.\n');

%% ===== Build pairs and calculate NORMALIZED coupling =====
pairs  = nchoosek(1:nClusters, 2);
nPairs = size(pairs, 1);
depthDiff = zeros(nPairs, 1);

% Zero-lag window in the SHORT CCG (for reference)
zeroBin_short = ceil(length(t)/2);
win_short     = zeroBin_short-2 : zeroBin_short+2;

% Zero-lag window in the LONG CCG (for normalized coupling)
zeroBin_long = ceil(length(t_long)/2);
win_long     = zeroBin_long-2 : zeroBin_long+2;

% Baseline bins: use bins beyond 100ms from zero lag in the long CCG
baselineBins_long = find(abs(t_long) >= 0.5); % 500ms

coupling_wake = zeros(nPairs, 1);
coupling_nrem = zeros(nPairs, 1);
coupling_rem  = zeros(nPairs, 1);

%% Calculate firing rates per unit per state (for filtering)
minFiringRate = 0.5; % Hz - minimum firing rate to be included

% Calculate duration of each state
wakeDuration = sum(wakeInts(:,2) - wakeInts(:,1));
nremDuration = sum(nremInts(:,2) - nremInts(:,1));
remDuration  = sum(remInts(:,2)  - remInts(:,1));

fprintf('State durations - WAKE: %.1f s, NREM: %.1f s, REM: %.1f s\n', ...
    wakeDuration, nremDuration, remDuration);

% Calculate firing rates
firingRate_wake = zeros(nClusters, 1);
firingRate_nrem = zeros(nClusters, 1);
firingRate_rem  = zeros(nClusters, 1);

for k = 1:nClusters
    firingRate_wake(k) = sum(s_wake(:,2) == k) / wakeDuration;
    firingRate_nrem(k) = sum(s_nrem(:,2) == k) / nremDuration;
    firingRate_rem(k)  = sum(s_rem(:,2)  == k) / remDuration;
end

fprintf('Units below %.1f Hz during WAKE: %d/%d\n', minFiringRate, sum(firingRate_wake < minFiringRate), nClusters);
fprintf('Units below %.1f Hz during NREM: %d/%d\n', minFiringRate, sum(firingRate_nrem < minFiringRate), nClusters);
fprintf('Units below %.1f Hz during REM:  %d/%d\n', minFiringRate, sum(firingRate_rem  < minFiringRate), nClusters);

%%

for p = 1:nPairs
    i = pairs(p,1);
    j = pairs(p,2);
    
    depthDiff(p) = abs(unitDepth(i) - unitDepth(j));
    
    % --- WAKE ---
    if firingRate_wake(i) < minFiringRate || firingRate_wake(j) < minFiringRate
        coupling_wake(p) = 0;
    else
        ccg_short_pair  = squeeze(ccg_wake(:,i,j));
        ccg_long_pair   = squeeze(ccg_wake_long(:,i,j));
        baseline        = median(ccg_long_pair(baselineBins_long));
        baseline_scaled = baseline * (binSize_short / binSize_long);
        if baseline_scaled == 0
            coupling_wake(p) = 0;
        else
            coupling_wake(p) = sum(ccg_short_pair(win_short)) / baseline_scaled - length(win_short);
        end
    end
    
    % --- NREM ---
    if firingRate_nrem(i) < minFiringRate || firingRate_nrem(j) < minFiringRate
        coupling_nrem(p) = 0;
    else
        ccg_short_pair  = squeeze(ccg_nrem(:,i,j));
        ccg_long_pair   = squeeze(ccg_nrem_long(:,i,j));
        baseline        = median(ccg_long_pair(baselineBins_long));
        baseline_scaled = baseline * (binSize_short / binSize_long);
        if baseline_scaled == 0
            coupling_nrem(p) = 0;
        else
            coupling_nrem(p) = sum(ccg_short_pair(win_short)) / baseline_scaled - length(win_short);
        end
    end
    
    % --- REM ---
    if firingRate_rem(i) < minFiringRate || firingRate_rem(j) < minFiringRate
        coupling_rem(p) = 0;
    else
        ccg_short_pair  = squeeze(ccg_rem(:,i,j));
        ccg_long_pair   = squeeze(ccg_rem_long(:,i,j));
        baseline        = median(ccg_long_pair(baselineBins_long));
        baseline_scaled = baseline * (binSize_short / binSize_long);
        if baseline_scaled == 0
            coupling_rem(p) = 0;
        else
            coupling_rem(p) = sum(ccg_short_pair(win_short)) / baseline_scaled - length(win_short);
        end
    end
end

fprintf('Normalized coupling done.\n');

%% Heatmap using LONG CCG
states   = {'WAKE', 'NREM', 'REM'};
all_ccgs_long = {ccg_wake_long, ccg_nrem_long, ccg_rem_long};
nStates  = length(states);
nUnits   = size(ccg_wake_long, 2);
nBins_long = size(ccg_wake_long, 1);

% Gaussian smoothing kernel
sigma = 3; % bins
xGauss = -3*sigma:3*sigma;
gaussKernel = exp(-xGauss.^2/(2*sigma^2));
gaussKernel = gaussKernel / sum(gaussKernel);

% Valid baseline bins for heatmap normalization
baselineWindow_heatmap = find(abs(t_long) >= 0.5); % beyond 500ms

%% Build pairOrder to track normCCG row ordering
pairOrder = zeros(nUnits*(nUnits-1)/2, 2);
rowCount  = 0;
for i = 1:nUnits
    for j = 1:i-1
        rowCount = rowCount + 1;
        pairOrder(rowCount, :) = [i j];
    end
end

%% Build normalized CCG matrix for heatmap
normCCG_states = cell(1, nStates);

for s = 1:nStates
    ccg    = all_ccgs_long{s};
    normCCG = [];
    
    for i = 1:nUnits
        for j = 1:i-1
            pairCCG = squeeze(ccg(:,i,j))';
            
            % Smooth
            pairCCG_sm = conv(pairCCG, gaussKernel, 'same');
            
            % Normalize using valid far-away baseline bins
            baseline = median(pairCCG_sm(baselineWindow_heatmap));
            
            if baseline == 0
                pairCCGnorm = pairCCG_sm;
            else
                pairCCGnorm = pairCCG_sm / baseline;
            end
            
            normCCG = [normCCG; pairCCGnorm];
        end
    end
    
    normCCG_states{s} = normCCG;
end

%% Map coupling_wake to normCCG row order for sorting
nNormRows = size(pairOrder, 1);
coupling_wake_mapped = zeros(nNormRows, 1);

for row = 1:nNormRows
    i = pairOrder(row, 1);
    j = pairOrder(row, 2);
    
    pairIdx = find( ...
        (pairs(:,1) == i & pairs(:,2) == j) | ...
        (pairs(:,1) == j & pairs(:,2) == i));
    
    if ~isempty(pairIdx)
        coupling_wake_mapped(row) = coupling_wake(pairIdx);
    end
end

%% Sort and plot heatmap
[~, sortIdx] = sort(coupling_wake_mapped, 'descend');

pairsToPlot = 400;
subsetIdx   = sortIdx(1:pairsToPlot);

figure;
for s = 1:nStates
    subplot(1, nStates, s)
    dataSubset = normCCG_states{s}(subsetIdx, :);
    imagesc(t_long, 1:pairsToPlot, dataSubset);
    xlabel('Time lag (s)');
    ylabel('Neuron pairs (sorted by WAKE coupling)');
    title([states{s} ' normalized CCGs']);
    clim([0.5 2.5]);
    colormap('hot');
    colorbar;
    set(gca, 'YDir', 'reverse');
    xlim([t_long(1) t_long(end)]);
end

% %% ===== HEATMAP (using short CCGs) =====
% states    = {'WAKE', 'NREM', 'REM'};
% all_ccgs  = {ccg_wake, ccg_nrem, ccg_rem};
% nStates   = length(states);
% nUnits    = size(ccg_wake, 2);
% nBins     = size(ccg_wake, 1);
% 
% % Gaussian smoothing kernel
% sigma = 3;
% xGauss = -3*sigma:3*sigma;
% gaussKernel = exp(-xGauss.^2/(2*sigma^2));
% gaussKernel = gaussKernel / sum(gaussKernel);
% 
% % Baseline window for heatmap normalization (use edge bins of short CCG)
% % Since short CCG is only 40ms, use outermost 20% of bins as baseline
% nEdgeBins = floor(0.2 * nBins);
% 
% normCCG_states = cell(1, nStates);
% 
% pairOrder = zeros(nUnits*(nUnits-1)/2, 2);
% rowCount = 0;
% for i = 1:nUnits
%     for j = 1:i-1
%         rowCount = rowCount + 1;
%         pairOrder(rowCount, :) = [i j]; % [unitA unitB]
%     end
% end
% 
% for s = 1:nStates
%     ccg = all_ccgs{s};
%     normCCG = [];
% 
%     for i = 1:nUnits
%         for j = 1:i-1
%             pairCCG = squeeze(ccg(:,i,j))';
% 
%             % Smooth
%             pairCCG_sm = conv(pairCCG, gaussKernel, 'same');
% 
%             % Normalize using edge bins of short CCG
%             baselineBins = [1:nEdgeBins, nBins-nEdgeBins+1:nBins];
%             baseline = median(pairCCG_sm(baselineBins));
% 
%             if baseline == 0
%                 pairCCGnorm = pairCCG_sm;
%             else
%                 pairCCGnorm = pairCCG_sm / baseline;
%             end
% 
%             normCCG = [normCCG; pairCCGnorm];
%         end
%     end
% 
%     normCCG_states{s} = normCCG;
% end
% 
% %% Map coupling_wake values to normCCG row order
% % For each row in normCCG, find the matching pair in the pairs matrix
% % and get its coupling value
% nNormRows = size(pairOrder, 1);
% coupling_wake_mapped = zeros(nNormRows, 1);
% 
% for row = 1:nNormRows
%     i = pairOrder(row, 1);
%     j = pairOrder(row, 2);
% 
%     % Find this pair in the pairs matrix (order might be [i,j] or [j,i])
%     pairIdx = find( ...
%         (pairs(:,1) == i & pairs(:,2) == j) | ...
%         (pairs(:,1) == j & pairs(:,2) == i) );
% 
%     if ~isempty(pairIdx)
%         coupling_wake_mapped(row) = coupling_wake(pairIdx);
%     end
% end
% 
% %% Sort by mapped coupling values
% [~, sortIdx] = sort(coupling_wake_mapped, 'descend');
% 
% pairsToPlot = length(sortIdx); %400;
% subsetIdx   = sortIdx(1:pairsToPlot);
% 
% 
% % % Build pair index mapping to match normCCG row ordering
% % % normCCG rows correspond to lower triangle pairs in order (i>j)
% % pairMap = zeros(nUnits*(nUnits-1)/2, 2);
% % rowCount = 0;
% % for i = 1:nUnits
% %     for j = 1:i-1
% %         rowCount = rowCount + 1;
% %         pairMap(rowCount,:) = [i j];
% %     end
% % end
% 
% %% Plot heatmaps
% figure;
% for s = 1:nStates
%     subplot(1, nStates, s)
%     dataSubset = normCCG_states{s}(subsetIdx,:);
%     imagesc(t, 1:pairsToPlot, dataSubset);
%     xlabel('Time lag (s)');
%     ylabel('Neuron pairs (sorted by WAKE coupling)');
%     title([states{s} ' normalized CCGs']);
%     clim([0.5 2.5]);
%     colormap('hot');
%     colorbar;
%     set(gca, 'YDir', 'reverse');
%     xlim([t(1) t(end)]);
% end

%% Plot CCG for one specific pair across all states
all_ccgs      = {ccg_wake,      ccg_nrem,      ccg_rem};
all_ccgs_long = {ccg_wake_long, ccg_nrem_long, ccg_rem_long};

clusterA = 89;%63
clusterB = 267; %415

unitA = find(good_clusters == clusterA);
unitB = find(good_clusters == clusterB);

if isempty(unitA) || isempty(unitB)
    error('One or both cluster IDs not found in good_clusters');
end

disp(['Cluster ' num2str(clusterA) ' → Unit ' num2str(unitA)])
disp(['Cluster ' num2str(clusterB) ' → Unit ' num2str(unitB)])

figure; hold on;
colors    = {'k','r','b'};
nBins    = size(ccg_wake, 1);
nEdgeBins = floor(0.2 * nBins); % edge bins of the short CCG

for s = 1:length(states)
    ccgPair         = squeeze(all_ccgs{s}(:, unitA, unitB));
    ccgLong         = squeeze(all_ccgs_long{s}(:, unitA, unitB));
    baseline        = median(ccgLong(baselineBins_long));
    baseline_scaled = baseline * (binSize_short / binSize_long);
    
    if baseline_scaled == 0
        ccgNorm = ccgPair;
    else
        ccgNorm = ccgPair / baseline_scaled;
    end
    
    plot(t, ccgNorm, 'Color', colors{s}, 'LineWidth', 2);
end

xlabel('Time lag (s)');
ylabel('Normalized CCG');
title(sprintf('Cluster %d → Cluster %d (Units %d → %d)', clusterA, clusterB, unitA, unitB));
legend(states, 'Location', 'best');
grid on;
xlim([t(1) t(end)]);

%% Save
ccgSavePath = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\NEW_ccgs-all-units-40ms-Lag';
save(fullfile(ccgSavePath, 'ccgs_all.mat'), ...
    'ccg_wake', 'ccg_nrem', 'ccg_rem', 't', ...
    'ccg_wake_long', 'ccg_nrem_long', 'ccg_rem_long', 't_long', ...
    'coupling_wake', 'coupling_nrem', 'coupling_rem', ...
    'pairs', 'depthDiff', 'good_clusters', 'unitDepth', ...
    'baselineBins_long', 'binSize_short', 'binSize_long', ...
    'nSpikes_wake', 'nSpikes_nrem', 'nSpikes_rem', 'minSpikes');

disp('CCGs and normalized coupling saved!');