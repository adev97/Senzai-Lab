% make_and_rank_ccgs_and_summary_figs.m

% new script that makes ccgs over 1s window with 1ms bins in order to
% identify which units fire together. Will use this to rank the most
% correlated units and display them using the short time windows

% this script works, can simply click run on the entire script and it will
% make and save the summary figs in your specified folder at the end!

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master")) 
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\analysis"))

sr = 30000;
nchan_probe = 384;
dtype = 'int16';

% kilosort directory
% ksDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Elissa_Belluccini\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';
ksDir = 'D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';

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

%% Calculate firing rates
wakeDuration = sum(wakeInts(:,2) - wakeInts(:,1));
nremDuration = sum(nremInts(:,2) - nremInts(:,1));
remDuration  = sum(remInts(:,2)  - remInts(:,1));

% fprintf('State durations (s) - WAKE: %.1f | NREM: %.1f | REM: %.1f\n', ...
%     wakeDuration, nremDuration, remDuration);

% Preallocate
minFiringRate = 0.5; % Hz

firingRate_wake = zeros(nClusters,1);
firingRate_nrem = zeros(nClusters,1);
firingRate_rem  = zeros(nClusters,1);

% Compute firing rate for each unit
for k = 1:nClusters
    
    if wakeDuration > 0
        firingRate_wake(k) = sum(s_wake(:,2) == k) / wakeDuration;
    end
    
    if nremDuration > 0
        firingRate_nrem(k) = sum(s_nrem(:,2) == k) / nremDuration;
    end
    
    if remDuration > 0
        firingRate_rem(k) = sum(s_rem(:,2) == k) / remDuration;
    end
end

fprintf('Units below %.2f Hz:\n', minFiringRate);
fprintf('  WAKE: %d / %d\n', sum(firingRate_wake < minFiringRate), nClusters);
fprintf('  NREM: %d / %d\n', sum(firingRate_nrem < minFiringRate), nClusters);
fprintf('  REM : %d / %d\n', sum(firingRate_rem  < minFiringRate), nClusters);

%% making CCGs for +- 1s
% use simply median normalization, still remove the below 0.5hz units and
% fill those pairs them with NaNs

binSize  = 0.001; % 1ms bins
duration = 2; % 2s total window

fprintf('Computing 1ms/2s CCGs...\n');
[ccg_wake, t] = CCG(s_wake(:,1), s_wake(:,2), 'binSize', binSize, 'duration', duration);
[ccg_nrem, ~] = CCG(s_nrem(:,1), s_nrem(:,2), 'binSize', binSize, 'duration', duration);
[ccg_rem,  ~] = CCG(s_rem(:,1),  s_rem(:,2),  'binSize', binSize, 'duration', duration);
fprintf('CCGs done.\n');

%% Save CCGs and t
ccgSavePath = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\ccgs-1ms-2s";
save(fullfile(ccgSavePath, 'ccgs.mat'), "ccg_wake", "ccg_nrem", "ccg_rem", "t");

load(fullfile(ccgSavePath, "ccgs.mat"))

%% Build pairs and calculate normalized CCGs
pairs  = nchoosek(1:nClusters, 2);
nPairs = size(pairs,1);

coupling_wake = nan(nPairs,1);
coupling_nrem = nan(nPairs,1);
coupling_rem  = nan(nPairs,1);

% Define windows
% zero lag window (±5 ms)
zeroBin = ceil(length(t)/2);
centerWin = zeroBin-5 : zeroBin+5;

% baseline = last 500 ms
baselineMask = abs(t) >= 0.5;

% Normalized Coupling
for p = 1:nPairs
    
    i = pairs(p,1);
    j = pairs(p,2);
    
    % ---- WAKE ----
    if firingRate_wake(i) >= minFiringRate && firingRate_wake(j) >= minFiringRate
        ccg_pair = squeeze(ccg_wake(:,i,j));
        baseline = median(ccg_pair(baselineMask));
        if baseline > 0
            centerVal = mean(ccg_pair(centerWin));
            coupling_wake(p) = centerVal / baseline - 1;
        else
            coupling_wake(p) = 0;
        end
    end
    
    % ---- NREM ----
    if firingRate_nrem(i) >= minFiringRate && firingRate_nrem(j) >= minFiringRate
        ccg_pair = squeeze(ccg_nrem(:,i,j));
        baseline = median(ccg_pair(baselineMask));
        if baseline > 0
            centerVal = mean(ccg_pair(centerWin));
            coupling_nrem(p) = centerVal / baseline - 1;
        else
            coupling_nrem(p) = 0;
        end
    end
    
    % ---- REM ----
    if firingRate_rem(i) >= minFiringRate && firingRate_rem(j) >= minFiringRate
        ccg_pair = squeeze(ccg_rem(:,i,j));
        baseline = median(ccg_pair(baselineMask));
        if baseline > 0
            centerVal = mean(ccg_pair(centerWin));
            coupling_rem(p) = centerVal / baseline - 1;
        else
            coupling_rem(p) = 0;
        end
    end
end

fprintf('Baseline-median normalized coupling complete.\n');


%% Rank pairs (most to least coupled based on wake)

% Only keep non-NaN pairs
validIdx = ~isnan(coupling_wake);
validPairs = pairs(validIdx, :);
validCoupling = coupling_wake(validIdx);

% Sort descending
[sortedCoupling, sortIdx] = sort(validCoupling, 'descend');
rankedPairs = validPairs(sortIdx, :);

% Display top 30
topN = min(30, length(sortedCoupling));
fprintf('Top %d most coupled pairs (unit i, unit j) with coupling value:\n', topN);
for k = 1:topN
    fprintf('Pair %d: (%d, %d) -> %.2f\n', k, rankedPairs(k,1), rankedPairs(k,2), sortedCoupling(k));
end

%% Make summary plots for the 30 highest coupled pairs)

% load needed files (meanWav and cell id type)
pngOutputDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\TEST-NEW_PAIR_SUMMARIES_COUPLED_UNITS_WITH_RF';
if ~exist(pngOutputDir,'dir'), mkdir(pngOutputDir); end

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

% Make summary plots for top 30 pairs during wake
% all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};   % your already normalized CCGs

% Reorder coupling vectors to match rankedPairs
coupling_wake_sorted = sortedCoupling;

coupling_nrem_sorted = coupling_nrem(validIdx);
coupling_nrem_sorted = coupling_nrem_sorted(sortIdx);

coupling_rem_sorted  = coupling_rem(validIdx);
coupling_rem_sorted  = coupling_rem_sorted(sortIdx);

all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};

sortBy = 'wake';
topN   = 30;

%% Load RF map data for this mouse
load("\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_RFMapping\RF_maps\Spike_Data\Mouse08_20251007_810to2250_RFmap_SpikeRate.mat")


%% Make summary figures for first 30 highest coupled pairs
makePairSummaryPNGsWORKS(all_ccgs, t, rankedPairs, good_clusters, ...
                         coupling_wake_sorted, coupling_nrem_sorted, coupling_rem_sorted, ...
                         unitDepth, meanWav, cell_type, ...
                         xpos, ypos, sr, ...
                         sortBy, topN, pngOutputDir, RFmap);


%% Comparision of receptive fields for ON and OFF separately between coupled pairs

%% 1. Compare full receptive field matrix, pearson correlation

topPairs = rankedPairs(1:topN,:);

RFsim_ON  = nan(topN,1);
RFsim_OFF = nan(topN,1);

for k = 1:topN
    
    i = topPairs(k,1);
    j = topPairs(k,2);
    
    RF_i_ON  = RFmap(i).ON;
    RF_j_ON  = RFmap(j).ON;
    
    RF_i_OFF = RFmap(i).OFF;
    RF_j_OFF = RFmap(j).OFF;
    
    % Remove NaNs if present
    valid_ON  = ~isnan(RF_i_ON)  & ~isnan(RF_j_ON);
    valid_OFF = ~isnan(RF_i_OFF) & ~isnan(RF_j_OFF);
    
    RFsim_ON(k)  = corr(RF_i_ON(valid_ON),  RF_j_ON(valid_ON));
    RFsim_OFF(k) = corr(RF_i_OFF(valid_OFF), RF_j_OFF(valid_OFF));
    
end




%% 2. Compare only distance between RF centers (peak positons)

%% 3. Compare RF overlap































% %% REDOING THIS PART
% %% Make short CCGs (0.04s total, 0.001s bins)
% 
% % parameters
% binSize_syn  = 0.001; % 1ms bins
% duration_syn = 0.040; % 8s total window
% 
% baselineThresh = 3; % preselection threshold
% jitterThresh = 5; % final synaptic threshold
% 
% jitterWindow = 0.01; % s
% nShuffle = 100;        % reduced for speed (originally 200)
% 
% fprintf('Computing 1ms/0.04s CCGs...\n');
% [ccg_wake, t_syn] = CCG(s_wake(:,1), s_wake(:,2), 'binSize', binSize_syn, 'duration', duration_syn);
% [ccg_nrem, ~] = CCG(s_nrem(:,1), s_nrem(:,2), 'binSize', binSize_syn, 'duration', duration_syn);
% [ccg_rem,  ~] = CCG(s_rem(:,1),  s_rem(:,2),  'binSize', binSize_syn, 'duration', duration_syn);
% fprintf('CCGs done.\n');
% 
% 
% %% define windows (for monosynaptic latencies)
% zeroBin = ceil(length(t_syn)/2);
% 
% postWin = zeroBin+1 : zeroBin+3;   % +1 to +3 ms
% preWin  = zeroBin-3 : zeroBin-1;   % -3 to -1 ms
% 
% baselineBins = find(abs(t_syn) >= 0.01 & abs(t_syn) <= 0.018);
% 
% %% Baseline Z score for all pairs
% 
% pairs  = nchoosek(1:nClusters, 2);
% nPairs = size(pairs, 1);
% depthDiff = zeros(nPairs, 1);
% 
% zeroIdx = find(t_syn == 0);  % index of 0 lag
% if isempty(zeroIdx)
%     [~, zeroIdx] = min(abs(t_syn)); % safety if exact 0 not found
% end
% 
% nPairs = size(pairs,1);
% 
% zBase_wake = nan(nPairs,1);
% zBase_nrem = nan(nPairs,1);
% zBase_rem  = nan(nPairs,1);
% 
% for p = 1:nPairs
% 
%     i = pairs(p,1);
%     j = pairs(p,2);
% 
%     % ================= WAKE =================
%     if firingRate_wake(i) >= minFiringRate && firingRate_wake(j) >= minFiringRate
% 
%         ccg_pair = squeeze(ccg_wake(:,i,j));
% 
%         mu = mean(ccg_pair(baselineBins));
%         sd = std(ccg_pair(baselineBins));
% 
%         if sd > 0
%             postPeak = sum(ccg_pair(postWin));
%             zBase_wake(p) = (postPeak - mu*length(postWin)) ...
%                 / (sd * sqrt(length(postWin)));
%         end
%     end
% 
%     % ================= NREM =================
%     if firingRate_nrem(i) >= minFiringRate && firingRate_nrem(j) >= minFiringRate
% 
%         ccg_pair = squeeze(ccg_nrem(:,i,j));
% 
%         mu = mean(ccg_pair(baselineBins));
%         sd = std(ccg_pair(baselineBins));
% 
%         if sd > 0
%             postPeak = sum(ccg_pair(postWin));
%             zBase_nrem(p) = (postPeak - mu*length(postWin)) ...
%                 / (sd * sqrt(length(postWin)));
%         end
%     end
% 
%     % ================= REM =================
%     if firingRate_rem(i) >= minFiringRate && firingRate_rem(j) >= minFiringRate
% 
%         ccg_pair = squeeze(ccg_rem(:,i,j));
% 
%         mu = mean(ccg_pair(baselineBins));
%         sd = std(ccg_pair(baselineBins));
% 
%         if sd > 0
%             postPeak = sum(ccg_pair(postWin));
%             zBase_rem(p) = (postPeak - mu*length(postWin)) ...
%                 / (sd * sqrt(length(postWin)));
%         end
%     end
% end
% 
% %% Select pairs that have a z score > 3
% cand_wake = find(zBase_wake > baselineThresh);
% cand_nrem = find(zBase_nrem > baselineThresh);
% cand_rem  = find(zBase_rem  > baselineThresh);
% 
% % %% Jitter function
% % jitter_spikes = @(ts) ts + (rand(size(ts)) - 0.5) * jitterWindow;
% 
% %% Jitter Correction for WAKE
% if isempty(gcp('nocreate'))
%     parpool;
% end
% 
% % takes around 10 min
% 
% fprintf('Computing jitter-corrected Z for WAKE, NREM, REM...\n');
% 
% % WAKE
% zJitter_wake = computeJitterZ(s_wake, cand_wake, pairs, postWin, binSize_syn, nShuffle, jitterWindow);
% % NREM
% zJitter_nrem = computeJitterZ(s_nrem, cand_nrem, pairs, postWin, binSize_syn, nShuffle, jitterWindow);
% % REM
% zJitter_rem  = computeJitterZ(s_rem,  cand_rem, pairs, postWin, binSize_syn, nShuffle, jitterWindow);
% 
% fprintf('All states complete.\n');
% 
% 
% %% Final synaptic pairs
% syn_wake = find(zJitter_wake > jitterThresh);
% syn_nrem = find(zJitter_nrem > jitterThresh);
% syn_rem  = find(zJitter_rem  > jitterThresh);
% 
% fprintf('Synaptic pairs:\nWAKE: %d\nNREM: %d\nREM: %d\n', ...
%     length(syn_wake), length(syn_nrem), length(syn_rem));
% 
% %% Depth differences for synaptic pairs
% depth_wake = abs(unitDepth(pairs(syn_wake,1)) - unitDepth(pairs(syn_wake,2)));
% depth_nrem = abs(unitDepth(pairs(syn_nrem,1)) - unitDepth(pairs(syn_nrem,2)));
% depth_rem  = abs(unitDepth(pairs(syn_rem,1)) - unitDepth(pairs(syn_rem,2)));
% 
% %% Overlap across states
% common_pairs = intersect(intersect(syn_wake, syn_nrem), syn_rem);
% fprintf('Pairs common across all states: %d\n', length(common_pairs));
% 
% %% Map syn_wake ids to cluster ids
% 
% synapticPairs_wake = pairs(syn_wake, :);
% clusterIDs_wake = good_clusters(synapticPairs_wake);
% 
% % clusterIDs contain the cluster numbers for units that have high coupling
% 
% %% Build summary pngs for all pairs that pass jitter test in wake
% mouseDir = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08";
% pntOutputDir = fullfile(mouseDir, 'synaptic_pair_plots');
% 
% % Load additional data needed for the function
% meanWaveformDir   = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters';
% classificationDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SC-RGC-classification-with-kmeans-all-channels-BEST';
% 
% load(fullfile(meanWaveformDir,   'meanWav_units.mat'));
% load(fullfile(classificationDir, 'RGC-SC-classification.mat'));
% 
% xpos = chanPos(:,1);
% ypos = chanPos(:,2);
% 
% % Flip cell_type orientation if needed
% cell_type  = cell_type';
% cluster_id = cluster_id';
% 
% 
% makePairSummaryPNGs(all_ccgs, all_ccgs_long, t_syn, pairs, syn_wake, good_clusters, ...
%                     zJitter_wake, zJitter_nrem, zJitter_rem, ...
%                     depthDiff, unitDepth, meanWav, cell_type, ...
%                     xpos, ypos, sr, baselineBins_long, binSize_short, binSize_long, ...
%                     pngDir);
% 















% %% Build pairs and calculate normalized coupling
% pairs  = nchoosek(1:nClusters, 2);
% nPairs = size(pairs, 1);
% depthDiff = zeros(nPairs, 1);
% 
% zeroIdx = find(t_syn == 0);  % index of 0 lag
% if isempty(zeroIdx)
%     [~, zeroIdx] = min(abs(t_syn)); % safety if exact 0 not found
% end
% 
% % preallocate
% coupling_wake = zeros(nPairs,1);
% coupling_nrem = zeros(nPairs,1);
% coupling_rem  = zeros(nPairs,1);
% 
% fprintf('Computing normalized coupling...\n');
% 
% for p = 1:nPairs
% 
%     i = pairs(p,1);
%     j = pairs(p,2);
% 
%     % raw zero-lag count
%     c_w = ccg_wake(zeroIdx,i,j);
%     c_n = ccg_nrem(zeroIdx,i,j);
%     c_r = ccg_rem(zeroIdx,i,j);
% 
%     % firing rates (for normalization)
%     fr_i_w = sum(s_wake(:,2)==i) / sum(diff(wakeInts,1,2));
%     fr_j_w = sum(s_wake(:,2)==j) / sum(diff(wakeInts,1,2));
% 
%     fr_i_n = sum(s_nrem(:,2)==i) / sum(diff(nremInts,1,2));
%     fr_j_n = sum(s_nrem(:,2)==j) / sum(diff(nremInts,1,2));
% 
%     fr_i_r = sum(s_rem(:,2)==i) / sum(diff(remInts,1,2));
%     fr_j_r = sum(s_rem(:,2)==j) / sum(diff(remInts,1,2));
% 
%     % expected coincidence under independence
%     exp_w = fr_i_w * fr_j_w * binSize_syn * sum(diff(wakeInts,1,2));
%     exp_n = fr_i_n * fr_j_n * binSize_syn * sum(diff(nremInts,1,2));
%     exp_r = fr_i_r * fr_j_r * binSize_syn * sum(diff(remInts,1,2));
% 
%     % normalized coupling (observed / expected)
%     if exp_w > 0
%         coupling_wake(p) = c_w / exp_w;
%     end
%     if exp_n > 0
%         coupling_nrem(p) = c_n / exp_n;
%     end
%     if exp_r > 0
%         coupling_rem(p) = c_r / exp_r;
%     end
% 
%     % depth difference
%     depthDiff(p) = abs(unitDepth(i) - unitDepth(j));
% 
% end
% 
% fprintf('Coupling computation done.\n');
% 
% 





