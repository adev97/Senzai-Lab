% FINAL_make_and_rank_ccgs_and_summary_figs.m

% new script that makes ccgs over 8s window with 1ms bins in order to
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

%% ------ filter low firing rate units out of subsequent calculations ---------
% Filter low-firing units (< 0.5 Hz in ANY state)

lowFiringUnits = (firingRate_wake < minFiringRate) | ...
                 (firingRate_nrem < minFiringRate) | ...
                 (firingRate_rem  < minFiringRate);

keepUnits = find(~lowFiringUnits);  % indices of units to keep
nClustersFiltered = length(keepUnits);

fprintf('Removing %d / %d units with FR < %.2f Hz in any state\n', ...
    sum(lowFiringUnits), nClusters, minFiringRate);

% Remap: filter good_clusters and unitDepth
% good_clusters = good_clusters(keepUnits);
% unitDepth     = unitDepth(keepUnits);

% Filter spike trains — remap unit IDs to new indices
keepMask = ismember(unitIDs, keepUnits);
ts = ts(keepMask);
oldIDs = unitIDs(keepMask);

% Remap old unit indices to new 1:nClustersFiltered indices
newIDs = zeros(size(oldIDs));
for i = 1:nClustersFiltered
    newIDs(oldIDs == keepUnits(i)) = i;
end

s = [ts newIDs];
nClusters = nClustersFiltered;

%% -------------------------------------------------------------------------- new

%% Load and verify external data alignment (before applying keepUnits filter)

% 1. Load and check meanWav
meanWaveformDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters';
good_clusters_main = good_clusters;          % save FIRST
load(fullfile(meanWaveformDir, 'meanWav_units.mat'));  % now load
good_clusters_wav  = good_clusters;          % rename loaded version
good_clusters      = good_clusters_main;     % restore original
% load(fullfile(meanWaveformDir, 'meanWav_units.mat'));
% good_clusters_main = good_clusters;          % save your main one before load overwrites it
% good_clusters_wav = good_clusters; % meanWav_units.mat saves good_clusters
% good_clusters      = good_clusters_main;     % restore your main one

% rename to avoid overwriting your main good_clusters
% check what variable name is saved in that mat file first

if ~isequal(good_clusters_wav, good_clusters)
    fprintf('WARNING: meanWav mismatch - reordering\n');
    % [~, reorderIdx] = ismember(good_clusters, good_clusters_wav);
    % if any(reorderIdx == 0); error('Units missing from meanWav'); end
    % meanWav = meanWav(:,:,reorderIdx);
else
    fprintf('meanWav OK\n');
end

% 2. Load and check cell_type
classificationDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SC-RGC-classification-with-kmeans-all-channels-BEST';
load(fullfile(classificationDir, 'RGC-SC-classification.mat'));
cell_type  = cell_type';
cluster_id = cluster_id';

if ~isequal(cluster_id, good_clusters)
    fprintf('WARNING: cell_type mismatch - reordering\n');
    % [~, reorderIdx] = ismember(good_clusters, cluster_id);
    % if any(reorderIdx == 0); error('Units missing from cell_type'); end
    % cell_type = cell_type(reorderIdx);
else
    fprintf('cell_type OK\n');
end

% 3. Load and check RFmap
load("\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_RFMapping\RF_maps\Spike_Data\Mouse08_20251007_810to2250_RFmap_SpikeRate.mat")
labels_file   = fullfile(ksDir, 'cluster_KSLabel.tsv');
unit_labels   = readtable(labels_file, 'FileType', 'text', 'Delimiter', '\t');
good_idx_rf   = strcmp(unit_labels.KSLabel, 'good');
good_units_rf = unit_labels.cluster_id(good_idx_rf);

if ~isequal(good_units_rf, good_clusters)
    fprintf('WARNING: RFmap mismatch - reordering\n');
    % [~, reorderIdx] = ismember(good_clusters, good_units_rf);
    % if any(reorderIdx == 0); error('Units missing from RFmap'); end
    % RFmap = RFmap(reorderIdx);
else
    fprintf('RFmap OK\n');
end

% 4. Now apply firing rate filter to all together
good_clusters = good_clusters(keepUnits);
unitDepth     = unitDepth(keepUnits);
meanWav       = meanWav(:,:,keepUnits);
cell_type     = cell_type(keepUnits);
RFmap         = RFmap(keepUnits);

%% --------------------------------------------------------------------------------------- new

% re-restrict the sleep states for the units used
% NOW re-restrict with the filtered spike train
s_wake = Restrict(s, wakeInts);
s_nrem = Restrict(s, nremInts);
s_rem  = Restrict(s, remInts);

%% making CCGs for +- 1s
% use simple median normalization, still remove the below 0.5hz units and
% fill those pairs them with NaNs

binSize  = 0.001; % 1ms bins
duration = 4; % 8s total window

tic

fprintf('Computing 1ms/8s CCGs...\n');
[ccg_wake, t] = CCG(s_wake(:,1), s_wake(:,2), 'binSize', binSize, 'duration', duration);
[ccg_nrem, ~] = CCG(s_nrem(:,1), s_nrem(:,2), 'binSize', binSize, 'duration', duration);
[ccg_rem,  ~] = CCG(s_rem(:,1),  s_rem(:,2),  'binSize', binSize, 'duration', duration);
fprintf('CCGs done.\n');

toc

all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};
states = {'WAKE', 'NREM', 'REM'};
nStates = length(states);

nBins = size(ccg_wake, 1);
[nTimeBins, nUnits, ~] = size(ccg_wake);

%% Save CCGs and t
ccgSavePath = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\passed-ccgs-1ms-8s";
if ~exist(ccgSavePath,'dir'), mkdir(ccgSavePath); end

save(fullfile(ccgSavePath, 'ccgs.mat'), "ccg_wake", "ccg_nrem", "ccg_rem", "t");

load(fullfile(ccgSavePath, "ccgs.mat"))

%% REDO THIS TO FIX NORMALIZATION Build pairs and calculate normalized CCGs
zeroBin = ceil(nTimeBins/2);

halfWidthBins = 2;   % ±2 bins = ±5 ms if binSize=1ms

centerWin = (zeroBin-halfWidthBins):(zeroBin+halfWidthBins);

pairs = nchoosek(1:nClusters, 2);
nPairs = size(pairs,1);

% depthDiff = zeros(nPairs,1);

all_ccgs_normalized = cell(size(all_ccgs));

coupling_wake = zeros(nUnits, nUnits);
coupling_nrem = zeros(nUnits, nUnits);
coupling_rem  = zeros(nUnits, nUnits);

% leftBins  = 1:(centerWin(1)-1);
% rightBins = (centerWin(end)+1):nTimeBins;
% baselineBins = [leftBins rightBins];
% 

for s = 1:nStates
    
    fprintf('Normalizing %s CCGs...\n', states{s});
    
    ccgState = all_ccgs{s};
    
    % allocate normalized matrix same size
    ccgNormState = zeros(size(ccgState));
    
    for pp = 1:nPairs
        
        a = pairs(pp,1);
        b = pairs(pp,2);
        
        % extract pair CCG
        ccgPair = squeeze(ccgState(:,a,b));
        
        % compute baseline from edges (excluding center)
        % baseline = median(ccgPair(baselineBins)); %%%
        baseline = median(ccgPair);   % FULL window median like yutas script

        % ----- Handle zero baseline properly -----
        if baseline == 0
            
            % No measurable coincidences → neutral baseline
            ccgNorm = ones(size(ccgPair));
            couplingValue = 1;
            
        else
            
            % Fold-change normalization
            ccgNorm = ccgPair / baseline;
            
            % Scalar coupling from center window
            couplingValue = mean(ccgNorm(centerWin));
            
        end
        
        % normalize FULL CCG
        % ccgNorm = ccgPair - baseline; % %%% subtraction baseline calculation (was division before and made huge blowups)
        % ccgNorm = ccgPair / baseline;

        % skip pairs with zero baseline to avoid blowups
        % if baseline == 0
        %     continue
        % end

        % store normalized CCG
        % ccgNorm = ccgPair / baseline;
        
        % store normalized CCG
        ccgNormState(:,a,b) = ccgNorm;
        ccgNormState(:,b,a) = flipud(ccgNorm); % maintain symmetry
        
        % compute scalar coupling from normalized CCG center window
        % couplingValue = mean(ccgNorm(centerWin));
        
        switch s
            case 1
                coupling_wake(a,b) = couplingValue;
                coupling_wake(b,a) = couplingValue;
            case 2
                coupling_nrem(a,b) = couplingValue;
                coupling_nrem(b,a) = couplingValue;
            case 3
                coupling_rem(a,b) = couplingValue;
                coupling_rem(b,a) = couplingValue;
        end
        
    end
    
    % store normalized state
    all_ccgs_normalized{s} = ccgNormState;
    
end

fprintf('Normalization and coupling calculation complete.\n');


%% ---------------- Comparision of whether units coupled in wake are coupled in nrem or rem
%% wake and nrem
validAll = ~isnan(coupling_wake) & ~isnan(coupling_nrem);

figure('Color', 'w', 'Position', [100, 100, 500, 500]);
scatter(coupling_wake(validAll), coupling_nrem(validAll), 20, 'filled', 'MarkerFaceAlpha', 0.4);
hold on;

% Add Unity Line (where Wake = NREM)
maxVal = max([coupling_wake(validAll); coupling_nrem(validAll)]);
minVal = min([coupling_wake(validAll); coupling_nrem(validAll)]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1.5);

xlabel('Wake Coupling Strength');
ylabel('NREM Coupling Strength');
title('Coupling Wake and NREM');
grid on; axis square;

% Calculate Correlation
[r, p] = corr(coupling_wake(validAll), coupling_nrem(validAll));
legend(sprintf('r = %.3f, p = %.3e', r, p), 'Unity Line', 'Location', 'northwest');

%% wake and rem
validAll = ~isnan(coupling_wake) & ~isnan(coupling_rem);

figure('Color', 'w', 'Position', [100, 100, 500, 500]);
scatter(coupling_wake(validAll), coupling_rem(validAll), 20, 'filled', 'MarkerFaceAlpha', 0.4);
hold on;

% Add Unity Line (where Wake = NREM)
maxVal = max([coupling_wake(validAll); coupling_rem(validAll)]);
minVal = min([coupling_wake(validAll); coupling_rem(validAll)]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1.5);

xlabel('Wake Coupling Strength');
ylabel('REM Coupling Strength');
title('Coupling Wake and REM');
grid on; axis square;

% Calculate Correlation
[r, p] = corr(coupling_wake(validAll), coupling_rem(validAll));
legend(sprintf('r = %.3f, p = %.3e', r, p), 'Unity Line', 'Location', 'northwest');

%% nrem and rem
validAll = ~isnan(coupling_rem) & ~isnan(coupling_nrem);

figure('Color', 'w', 'Position', [100, 100, 500, 500]);
scatter(coupling_nrem(validAll), coupling_rem(validAll), 20, 'filled', 'MarkerFaceAlpha', 0.4);
hold on;

% Add Unity Line (where Wake = NREM)
maxVal = max([coupling_nrem(validAll); coupling_rem(validAll)]);
minVal = min([coupling_nrem(validAll); coupling_rem(validAll)]);
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1.5);

xlabel('NREM Coupling Strength');
ylabel('REM Coupling Strength');
title('Coupling NREM and REM');
grid on; axis square;

% Calculate Correlation
[r, p] = corr(coupling_nrem(validAll), coupling_rem(validAll));
legend(sprintf('r = %.3f, p = %.3e', r, p), 'Unity Line', 'Location', 'northwest');


%% ----------- Rank pairs (most to least coupled based on wake)

nPairs = size(pairs,1);

% Extract coupling values corresponding to each pair
couplingWakeVec = nan(nPairs,1);
couplingNremVec = nan(nPairs,1);
couplingRemVec  = nan(nPairs,1);

for k = 1:nPairs
    a = pairs(k,1);
    b = pairs(k,2);

    couplingWakeVec(k) = coupling_wake(a,b);
    couplingNremVec(k) = coupling_nrem(a,b);
    couplingRemVec(k)  = coupling_rem(a,b);
end

% Only keep valid pairs
validIdx = ~isnan(couplingWakeVec);

validPairs       = pairs(validIdx,:);
validWake        = couplingWakeVec(validIdx);
validNrem        = couplingNremVec(validIdx);
validRem         = couplingRemVec(validIdx);

% Sort descending by wake coupling
[sortedWake, sortIdx] = sort(validWake, 'descend');

rankedPairs = validPairs(sortIdx,:);

coupling_wake_sorted = sortedWake;
coupling_nrem_sorted = validNrem(sortIdx);
coupling_rem_sorted  = validRem(sortIdx);

% Display top pairs
topN = min(30, length(sortedWake));

fprintf('Top %d most coupled pairs:\n', topN);

for k = 1:topN
    fprintf('Pair %d: (%d, %d) -> Wake=%.3f, NREM=%.3f, REM=%.3f\n', ...
        k, ...
        rankedPairs(k,1), ...
        rankedPairs(k,2), ...
        coupling_wake_sorted(k), ...
        coupling_nrem_sorted(k), ...
        coupling_rem_sorted(k));
end

%% ------- Make summary plots for the 30 highest coupled pairs)

% load needed files (meanWav and cell id type)
pngOutputDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\hmm-TEST3-8s-plot-1s_with_norm';
if ~exist(pngOutputDir,'dir'), mkdir(pngOutputDir); end

% Load additional data needed for the function
% meanWaveformDir   = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters';
% classificationDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SC-RGC-classification-with-kmeans-all-channels-BEST';
% 
% load(fullfile(meanWaveformDir,   'meanWav_units.mat'));
% load(fullfile(classificationDir, 'RGC-SC-classification.mat'));

% Flip cell_type orientation if needed
% cell_type  = cell_type';
% cluster_id = cluster_id';

% Make summary plots for top 30 pairs during wake
% all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};   % your already normalized CCGs

% Reorder coupling vectors to match rankedPairs
% coupling_wake_sorted = sortedCoupling;
% 
% coupling_nrem_sorted = coupling_nrem(validIdx);
% coupling_nrem_sorted = coupling_nrem_sorted(sortIdx);
% 
% coupling_rem_sorted  = coupling_rem(validIdx);
% coupling_rem_sorted  = coupling_rem_sorted(sortIdx);

% norm_all_ccgs = {coupling_wake, coupling_nrem, coupling_rem};

sortBy = 'wake';
topN   = 200;

%% Load RF map data for this mouse
% load("\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_RFMapping\RF_maps\Spike_Data\Mouse08_20251007_810to2250_RFmap_SpikeRate.mat")
% 
% if ~isequal(good_units, good_clusters)
%     fprintf('WARNING: good_units and good_clusters do not match - reordering RFmap\n');
% 
%     % Reorder RFmap to match good_clusters ordering
%     [~, reorderIdx] = ismember(good_clusters, good_units);
% 
%     if any(reorderIdx == 0)
%         error('Some units in good_clusters not found in good_units - RF mapping may be incomplete');
%     end
% 
%     RFmap = RFmap(reorderIdx);
% else
%     fprintf('good_units and good_clusters match - RFmap ordering is correct\n');
% end




% meanWav       = meanWav(:,:,keepUnits);
% cell_type     = cell_type(keepUnits);
% RFmap         = RFmap(keepUnits);

% sanity check
topCheck = 5;

figure('Color','w');
for k = 1:topCheck
    
    a = rankedPairs(k,1);
    b = rankedPairs(k,2);
    
    subplot(topCheck,1,k)
    
    ccgWake = squeeze(all_ccgs_normalized{1}(:,a,b)); % 1 = wake
    
    plot(t, ccgWake, 'k','LineWidth',1.5); hold on
    xline(0,'r--')
    
    title(sprintf('Rank %d | Units (%d,%d) | Wake coupling = %.3f', ...
        k, a, b, coupling_wake_sorted(k)))
    
    ylabel('Norm CCG')
    
    if k==topCheck
        xlabel('Time (ms)')
    end
    
end


%% Make summary figures for first 30 highest coupled pairs
makePairSummaryPNGsWORKS(all_ccgs_normalized, t, rankedPairs, good_clusters, ...
                         coupling_wake_sorted, coupling_nrem_sorted, coupling_rem_sorted, ...
                         unitDepth, meanWav, cell_type, ...
                         xpos, ypos, sr, ...
                         sortBy, topN, pngOutputDir, RFmap);


%% ------- Comparision of receptive fields for ON and OFF separately between coupled pairs

%% 1. Compare full receptive field matrix, pearson correlation -- results in a lot of zeros

topPairs = rankedPairs(1:topN,:);

rf_similarity_on  = zeros(topN,1);
rf_similarity_off = zeros(topN,1);

for k = 1:topN

    uA = topPairs(k,1);
    uB = topPairs(k,2);

    rfA_on  = mean(RFmap{uA}.ON.OnSet,3)  - RFmap{uA}.baseline;
    rfB_on  = mean(RFmap{uB}.ON.OnSet,3)  - RFmap{uB}.baseline;

    rfA_off = mean(RFmap{uA}.OFF.OnSet,3) - RFmap{uA}.baseline;
    rfB_off = mean(RFmap{uB}.OFF.OnSet,3) - RFmap{uB}.baseline;

    rf_similarity_on(k)  = corr(rfA_on(:),  rfB_on(:));
    rf_similarity_off(k) = corr(rfA_off(:), rfB_off(:));

end

figure;
scatter(rf_similarity_on, coupling_wake_sorted(1:topN))
xlabel('RF similarity (ON)')
ylabel('Coupling (wake)')

figure;
scatter(rf_similarity_off, coupling_wake_sorted(1:topN))
xlabel('RF similarity (OFF)')
ylabel('Coupling (wake)')

%% 2. Compare only distance between RF centers (peak positons)
% % use abs to make it polarity independent
% [~, idxA] = max(abs(rfA_on(:)));
% [yA, xA]  = ind2sub(size(rfA_on), idxA);
% 
% [~, idxB] = max(abs(rfB_on(:)));
% [yB, xB]  = ind2sub(size(rfB_on), idxB);
% 
% % compute peak distance (by pixel)
% peakDist = sqrt((xA - xB)^2 + (yA - yB)^2);
% 
% thresholdA = prctile(abs(rfA_on(:)), 95);
% maskA = abs(rfA_on) > thresholdA;
% 
% thresholdB = prctile(abs(rfB_on(:)), 95);
% maskB = abs(rfB_on) > thresholdB;
% 
% overlap = sum(maskA(:) & maskB(:)) / sum(maskA(:) | maskB(:));

%% 3. Compare RF overlap (Jaccard, simple Index)

% Parameters
nPairsToRun = size(rankedPairs, 1);
topPairs = rankedPairs(1:topN,:);
sigma = 1; % Smoothing factor for RF maps

% Preallocate
rf_overlap_on   = zeros(nPairsToRun,1);
rf_overlap_off  = zeros(nPairsToRun,1);
centroid_dist_on   = zeros(nPairsToRun,1); 
centroid_dist_off   = zeros(nPairsToRun,1);

maskVal = 90;

for k = 1:nPairsToRun
    uA = rankedPairs(k,1);
    uB = rankedPairs(k,2);
    
    % 1. Process RFs: Mean across time, subtract baseline, and Smooth
    % Smoothing helps the thresholding identify 'blobs' rather than single pixels
    rfA_on  = imgaussfilt(mean(RFmap{uA}.ON.OnSet,3)  - RFmap{uA}.baseline, sigma);
    rfB_on  = imgaussfilt(mean(RFmap{uB}.ON.OnSet,3)  - RFmap{uB}.baseline, sigma);
    rfA_off = imgaussfilt(mean(RFmap{uA}.OFF.OnSet,3) - RFmap{uA}.baseline, sigma);
    rfB_off = imgaussfilt(mean(RFmap{uB}.OFF.OnSet,3) - RFmap{uB}.baseline, sigma);
    
    % 2. Thresholding (Top 10% of signal)
    threshA_on = prctile(rfA_on(:), maskVal);
    threshB_on = prctile(rfB_on(:), maskVal);
    threshA_off = prctile(rfA_off(:), maskVal);
    threshB_off = prctile(rfB_off(:), maskVal);

    noiseThresh_A_on  = 2 * std(rfA_on(:)); % kep pixels that are two std from baseline firing
    noiseThresh_B_on  = 2 * std(rfB_on(:));
    noiseThresh_A_off = 2 * std(rfA_off(:));
    noiseThresh_B_off = 2 * std(rfB_off(:));

    maskA_on  = (rfA_on  >= threshA_on)  & (rfA_on  > noiseThresh_A_on);
    maskB_on  = (rfB_on  >= threshB_on)  & (rfB_on  > noiseThresh_B_on);
    maskA_off = (rfA_off >= threshA_off) & (rfA_off > noiseThresh_A_off);
    maskB_off = (rfB_off >= threshB_off) & (rfB_off > noiseThresh_B_off);
    
    % maskA_on = rfA_on >= threshA_on;
    % maskB_on = rfB_on >= threshB_on;
    % maskA_off = rfA_off >= threshA_off;
    % maskB_off = rfB_off >= threshB_off;

    % 3. Calculate Jaccard Overlap (Intersection over Union)
    intersection = sum(maskA_on(:) & maskB_on(:));
    union_area_on   = sum(maskA_on(:) | maskB_on(:));
    if union_area_on > 0
        rf_overlap_on(k) = intersection / union_area_on;
    else
        rf_overlap_on(k) = 0;
    end

    intersection_off = sum(maskA_off(:) & maskB_off(:));
    union_area_off   = sum(maskA_off(:) | maskB_off(:));
    if union_area_off > 0
        rf_overlap_off(k) = intersection_off / union_area_off;
    else
        rf_overlap_off(k) = 0;
    end
    
    % 4. Calculate Centroid Distance (The 'Center' of their representation)
    % We use regionprops on the binary mask to find the center of mass
    statsA_on = regionprops(maskA_on, 'Centroid');
    statsB_on = regionprops(maskB_on, 'Centroid');
    
    if ~isempty(statsA_on) && ~isempty(statsB_on)
        % Take the largest blob if multiple exist
        cA = statsA_on(1).Centroid; 
        cB = statsB_on(1).Centroid;
        centroid_dist_on(k) = sqrt((cA(1)-cB(1))^2 + (cA(2)-cB(2))^2);
    else
        centroid_dist_on(k) = NaN;
    end

    statsA_off = regionprops(maskA_off, 'Centroid');
    statsB_off = regionprops(maskB_off, 'Centroid');
    
    if ~isempty(statsA_off) && ~isempty(statsB_off)
        % Take the largest blob if multiple exist
        cA = statsA_off(1).Centroid; 
        cB = statsB_off(1).Centroid;
        centroid_dist_off(k) = sqrt((cA(1)-cB(1))^2 + (cA(2)-cB(2))^2);
    else
        centroid_dist_off(k) = NaN;
    end

end

%% 4. Compare RF Overlap (2D Gaussian)
% each receptive field gets one gaussian curve. then, compute the integral
% of the product of the two gaussians divided by integral of their union

% Preallocate Gaussian overlap
rf_overlap_on_gauss  = zeros(nPairsToRun,1);
rf_overlap_off_gauss = zeros(nPairsToRun,1);

tic
for k = 1:nPairsToRun
    uA = rankedPairs(k,1);
    uB = rankedPairs(k,2);

    % ----- ON RFs -----
    rfA_on = mean(RFmap{uA}.ON.OnSet,3) - RFmap{uA}.baseline;
    rfB_on = mean(RFmap{uB}.ON.OnSet,3) - RFmap{uB}.baseline;
    
    % Fit 2D Gaussian to each RF
    paramsA_on = fit2DGaussian(rfA_on);
    paramsB_on = fit2DGaussian(rfB_on);
    
    % Generate Gaussian surfaces
    [X,Y] = meshgrid(1:size(rfA_on,2), 1:size(rfA_on,1));
    GA = paramsA_on.A * exp(-(((X-paramsA_on.x0).^2)/(2*paramsA_on.sigma_x^2) + ((Y-paramsA_on.y0).^2)/(2*paramsA_on.sigma_y^2)));
    GB = paramsB_on.A * exp(-(((X-paramsB_on.x0).^2)/(2*paramsB_on.sigma_x^2) + ((Y-paramsB_on.y0).^2)/(2*paramsB_on.sigma_y^2)));
    
    % Overlap = integral(product) / integral(union)
    rf_overlap_on_gauss(k) = sum(GA(:).*GB(:)) / (sum(GA(:)) + sum(GB(:)) - sum(GA(:).*GB(:)));
    
    % ----- OFF RFs -----
    rfA_off = mean(RFmap{uA}.OFF.OnSet,3) - RFmap{uA}.baseline;
    rfB_off = mean(RFmap{uB}.OFF.OnSet,3) - RFmap{uB}.baseline;
    
    paramsA_off = fit2DGaussian(rfA_off);
    paramsB_off = fit2DGaussian(rfB_off);
    
    GA = paramsA_off.A * exp(-(((X-paramsA_off.x0).^2)/(2*paramsA_off.sigma_x^2) + ((Y-paramsA_off.y0).^2)/(2*paramsA_off.sigma_y^2)));
    GB = paramsB_off.A * exp(-(((X-paramsB_off.x0).^2)/(2*paramsB_off.sigma_x^2) + ((Y-paramsB_off.y0).^2)/(2*paramsB_off.sigma_y^2)));
    
    rf_overlap_off_gauss(k) = sum(GA(:).*GB(:)) / (sum(GA(:)) + sum(GB(:)) - sum(GA(:).*GB(:)));
end
toc

%% ----------- Make scatter for all units (rf overlap vs coupling strength)

% % Choose overlap metric
% useGaussianRF = true;  % true = Gaussian fit overlap, false = Jaccard index
% 
% if useGaussianRF
%     rf_valid    = rf_overlap_on_gauss(validIdx);   % Gaussian-based ON overlap
%     rf_valid_off= rf_overlap_off_gauss(validIdx);  % Gaussian-based OFF overlap
%     overlapLabel = 'Gaussian RF Overlap';
% else
%     rf_valid    = rf_overlap_on(validIdx);   % Original Jaccard ON
%     rf_valid_off= rf_overlap_off(validIdx);  % Original Jaccard OFF
%     overlapLabel = 'Jaccard RF Overlap';
% end
% 
% coupling_valid = coupling_wake_sorted(validIdx);
% 
% % Sort by coupling strength (descending)
% [sortedCoupling, sortIdx] = sort(coupling_valid, 'descend');
% rf_sorted     = rf_valid(sortIdx);
% rf_sorted_off = rf_valid_off(sortIdx);
% 
% num_to_plot = size(sortedCoupling, 1);  % plot all pairs
% x_data     = rf_sorted(1:num_to_plot);
% x_data_off = rf_sorted_off(1:num_to_plot);
% y_data     = sortedCoupling(1:num_to_plot);
% 
% %% ---- ON RF Overlap Scatter ----
% [r, p] = corr(x_data, y_data);
% 
% figure('Color','w','Position',[100,100,800,800]);
% scatter(x_data, y_data, 20, 'filled', 'MarkerFaceAlpha',0.3,'MarkerFaceColor',[0 0.45 0.74]);
% hold on;
% 
% % Trend line
% coeffs = polyfit(x_data, y_data, 1);
% fitX = linspace(min(x_data), max(x_data), 100);
% fitY = polyval(coeffs, fitX);
% plot(fitX, fitY, 'k', 'LineWidth', 1);
% 
% xlabel([overlapLabel ' (ON)']);
% ylabel('Wake Coupling Strength');
% title(sprintf('Functional Coupling vs. Spatial Overlap (ON)\nr = %.3f, p = %.2e', r, p));
% grid on;
% % xlim([-100 100]);
% 
% %% ---- OFF RF Overlap Scatter ----
% [r, p] = corr(x_data_off, y_data);
% 
% figure('Color','w','Position',[100,100,800,800]);
% scatter(x_data_off, y_data, 20, 'filled', 'MarkerFaceAlpha',0.3,'MarkerFaceColor',[0.85 0.33 0.1]);
% hold on;
% 
% % Trend line
% coeffs = polyfit(x_data_off, y_data, 1);
% fitX = linspace(min(x_data_off), max(x_data_off), 100);
% fitY = polyval(coeffs, fitX);
% plot(fitX, fitY, 'k', 'LineWidth', 1);
% 
% xlabel([overlapLabel ' (OFF)']);
% ylabel('Wake Coupling Strength');
% title(sprintf('Functional Coupling vs. Spatial Overlap (OFF)\nr = %.3f, p = %.2e', r, p));
% grid on;

% remove NaNs first
validIdx = ~isnan(rf_overlap_on) & ~isnan(coupling_wake_sorted);

rf_valid  = rf_overlap_on(validIdx);
rf_valid_off = rf_overlap_off(validIdx);
coupling_valid = coupling_wake_sorted(validIdx);

% sort by coupling strength (descending)
[sortedCoupling, sortIdx] = sort(coupling_valid, 'descend');

rf_sorted = rf_valid(sortIdx);
rf_sorted_off = rf_valid_off(sortIdx);

% % select top 50 coupled pairs
num_to_plot = min(500, length(sortedCoupling));

% if want to plot all
% num_to_plot = size(sortedCoupling, 1);

x_data     = rf_sorted(1:num_to_plot);
x_data_off = rf_sorted_off(1:num_to_plot);
y_data     = sortedCoupling(1:num_to_plot);

%% FOR ON
% Calculate Correlation
[r, p] = corr(x_data, y_data);

figure('Color', 'w', 'Position', [100, 100, 800, 800]);
% Use 'binscatter' if you have many points, or 'scatter' with transparency
scatter(x_data, y_data, 20, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0 0.45 0.74]);
hold on;

% Add a trend line (Linear Regression)
coeffs = polyfit(x_data, y_data, 1);
fitX = linspace(min(x_data), max(x_data), 100);
fitY = polyval(coeffs, fitX);
plot(fitX, fitY, 'k', 'LineWidth', 1);

xlabel('ON RF Overlap (Jaccard Index)');
ylabel('Wake Coupling Strength');
title(sprintf('Functional Coupling vs. Spatial Overlap\nr = %.3f (p = %.2e)', r, p));
grid on;

%% FOR OFF
% Calculate Correlation
[r, p] = corr(x_data_off, y_data);

figure('Color', 'w', 'Position', [100, 100, 800, 800]);
% Use 'binscatter' if you have many points, or 'scatter' with transparency
scatter(x_data_off, y_data, 20, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0.85 0.33 0.1]);
hold on;

% Add a trend line (Linear Regression)
coeffs = polyfit(x_data_off, y_data, 1);
fitX = linspace(min(x_data_off), max(x_data_off), 100);
fitY = polyval(coeffs, fitX);
plot(fitX, fitY, 'k', 'LineWidth', 1);

xlabel('OFF RF Overlap (Jaccard Index)');
ylabel('Wake Coupling Strength');
title(sprintf('Functional Coupling vs. Spatial Overlap\nr = %.3f (p = %.2e)', r, p));
grid on;

%% Test plot ON RF overlap vs OFF RF overlap
[r, p] = corr(x_data, x_data_off);
figure;
scatter(x_data, x_data_off, 20, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [1.0, 0.4, 0.6]);
hold on;
% Add a trend line (Linear Regression)
coeffs = polyfit(x_data, x_data_off, 1);
fitX = linspace(min(x_data), max(x_data), 100);
fitY = polyval(coeffs, fitX);
plot(fitX, fitY, 'k', 'LineWidth', 1);

xlabel('ON RF Overlap (Jaccard Index)');
ylabel('OFF RF Overlap (Jaccard Index)');
grid on;

title(sprintf('ON rf vs OFF rf\nr = %.3f (p = %.2e)', r, p));

%% ---- plot NREM coupling vs receptive field overlap ----
% NREM coupling vs RF overlap OFF
validIdx_nrem = ~isnan(rf_overlap_on) & ~isnan(coupling_nrem_sorted);

[r, p] = corr(rf_overlap_on(validIdx_nrem), coupling_nrem_sorted(validIdx_nrem));

figure('Color', 'w', 'Position', [100, 100, 800, 800]);
scatter(rf_overlap_on(validIdx_nrem), coupling_nrem_sorted(validIdx_nrem), 20, 'filled', ...
    'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0 0.45 0.74]);
hold on;
coeffs = polyfit(rf_overlap_on(validIdx_nrem), coupling_nrem_sorted(validIdx_nrem), 1);
fitX = linspace(min(rf_overlap_on(validIdx_nrem)), max(rf_overlap_on(validIdx_nrem)), 100);
plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1);
xlabel('ON RF Overlap (Jaccard Index)');
ylabel('NREM Coupling Strength');
title(sprintf('NREM Coupling vs. Spatial Overlap (ON)\nr = %.3f (p = %.2e)', r, p));
grid on;

% NREM coupling vs RF overlap OFF
validIdx_nrem = ~isnan(rf_overlap_off) & ~isnan(coupling_nrem_sorted);

[r, p] = corr(rf_overlap_off(validIdx_nrem), coupling_nrem_sorted(validIdx_nrem));

figure('Color', 'w', 'Position', [100, 100, 800, 800]);
scatter(rf_overlap_off(validIdx_nrem), coupling_nrem_sorted(validIdx_nrem), 20, 'filled', ...
    'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0.85 0.33 0.1]);
hold on;
coeffs = polyfit(rf_overlap_off(validIdx_nrem), coupling_nrem_sorted(validIdx_nrem), 1);
fitX = linspace(min(rf_overlap_off(validIdx_nrem)), max(rf_overlap_off(validIdx_nrem)), 100);
plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1);
xlabel('OFF RF Overlap (Jaccard Index)');
ylabel('NREM Coupling Strength');
title(sprintf('NREM Coupling vs. Spatial Overlap (OFF)\nr = %.3f (p = %.2e)', r, p));
grid on;

%% ----- plot REM coupling vs receptive field overlap ----
% REM coupling vs RF overlap ON
validIdx_rem = ~isnan(rf_overlap_on) & ~isnan(coupling_rem_sorted);

[r, p] = corr(rf_overlap_on(validIdx_rem), coupling_rem_sorted(validIdx_rem));

figure('Color', 'w', 'Position', [100, 100, 800, 800]);
scatter(rf_overlap_on(validIdx_rem), coupling_rem_sorted(validIdx_rem), 20, 'filled', ...
    'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0 0.45 0.74]);
hold on;
coeffs = polyfit(rf_overlap_on(validIdx_rem), coupling_rem_sorted(validIdx_rem), 1);
fitX = linspace(min(rf_overlap_on(validIdx_rem)), max(rf_overlap_on(validIdx_rem)), 100);
plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1);
xlabel('ON RF Overlap (Jaccard Index)');
ylabel('REM Coupling Strength');
title(sprintf('REM Coupling vs. Spatial Overlap (ON)\nr = %.3f (p = %.2e)', r, p));
grid on;

% REM coupling vs RF overlap OFF
validIdx_rem = ~isnan(rf_overlap_off) & ~isnan(coupling_rem_sorted);

[r, p] = corr(rf_overlap_off(validIdx_rem), coupling_rem_sorted(validIdx_rem));

figure('Color', 'w', 'Position', [100, 100, 800, 800]);
scatter(rf_overlap_off(validIdx_rem), coupling_rem_sorted(validIdx_rem), 20, 'filled', ...
    'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0.85 0.33 0.1]);
hold on;
coeffs = polyfit(rf_overlap_off(validIdx_rem), coupling_rem_sorted(validIdx_rem), 1);
fitX = linspace(min(rf_overlap_off(validIdx_rem)), max(rf_overlap_off(validIdx_rem)), 100);
plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1);
xlabel('OFF RF Overlap (Jaccard Index)');
ylabel('REM Coupling Strength');
title(sprintf('REM Coupling vs. Spatial Overlap (OFF)\nr = %.3f (p = %.2e)', r, p));
grid on;

%% Combined figure with all ON/OFF vs wake/rem/nrem
figure('Color', 'w', 'Position', [100, 100, 1400, 900]);

stateLabels   = {'WAKE', 'NREM', 'REM'};
couplingAll   = {coupling_wake_sorted, coupling_nrem_sorted, coupling_rem_sorted};
stateColors   = {[0 0 0], [0.8 0.2 0.2], [0.2 0.4 0.8]};

for col = 1:3  % columns = states (WAKE, NREM, REM)

    coupling = couplingAll{col};

    % --- ON row (row 1) ---
    subplot(2, 3, col);
    validIdx = ~isnan(rf_overlap_on) & ~isnan(coupling);
    x = rf_overlap_on(validIdx);
    y = coupling(validIdx);

    [r, p] = corr(x, y);
    scatter(x, y, 15, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', stateColors{col});
    hold on;
    coeffs = polyfit(x, y, 1);
    fitX = linspace(min(x), max(x), 100);
    plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1.2);
    xlabel('ON RF Overlap (Jaccard)');
    ylabel([stateLabels{col} ' Coupling']);
    title(sprintf('%s — ON\nr = %.3f, p = %.2e', stateLabels{col}, r, p));
    grid on;

    % --- OFF row (row 2) ---
    subplot(2, 3, col + 3);
    validIdx = ~isnan(rf_overlap_off) & ~isnan(coupling);
    x = rf_overlap_off(validIdx);
    y = coupling(validIdx);

    [r, p] = corr(x, y);
    scatter(x, y, 15, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', stateColors{col});
    hold on;
    coeffs = polyfit(x, y, 1);
    fitX = linspace(min(x), max(x), 100);
    plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1.2);
    xlabel('OFF RF Overlap (Jaccard)');
    ylabel([stateLabels{col} ' Coupling']);
    title(sprintf('%s — OFF\nr = %.3f, p = %.2e', stateLabels{col}, r, p));
    grid on;

end

sgtitle('RF Overlap vs Coupling Strength Across States', 'FontSize', 16, 'FontWeight', 'bold');

% Save
saveas(gcf, fullfile(pngOutputDir, 'RF_overlap_vs_coupling_all_states.png'));



%% New function call include overlap RFs in the summary pdf
pngOutputDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\0-PAIR-SUMMARIES-WITH-RF';
if ~exist(pngOutputDir,'dir'), mkdir(pngOutputDir); end

dispTime = 0.5; % amount of +- time plotted as the ccg in the summary figure

makePairSummaryPNGsPlusCorrRFs(all_ccgs_normalized, t, rankedPairs, good_clusters, ...
                         coupling_wake_sorted, coupling_nrem_sorted, coupling_rem_sorted, ...
                         unitDepth, meanWav, cell_type, ...
                         xpos, ypos, sr, ...
                         sortBy, topN, pngOutputDir, RFmap, maskVal, dispTime);









