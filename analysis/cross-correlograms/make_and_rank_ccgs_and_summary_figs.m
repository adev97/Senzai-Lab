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
% ccgSavePath = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\ccgs-1ms-2s";
% save(fullfile(ccgSavePath, 'ccgs.mat'), "ccg_wake", "ccg_nrem", "ccg_rem", "t");
% 
% load(fullfile(ccgSavePath, "ccgs.mat"))

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

%% ------- Make summary plots for the 30 highest coupled pairs)

% load needed files (meanWav and cell id type)
pngOutputDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SUM-NEW_PAIR_SUMMARIES_COUPLED_UNITS_WITH_RF';
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
topN   = 60;

%% Load RF map data for this mouse
load("\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_RFMapping\RF_maps\Spike_Data\Mouse08_20251007_810to2250_RFmap_SpikeRate.mat")


%% Make summary figures for first 30 highest coupled pairs
makePairSummaryPNGsWORKS(all_ccgs, t, rankedPairs, good_clusters, ...
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
    
    rf_similarity_on(i)  = corr(rfA_on(:),  rfB_on(:));
    rf_similarity_off(i) = corr(rfA_off(:), rfB_off(:));
    
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

%% 3. Compare RF overlap

% Parameters
nPairsToRun = size(rankedPairs, 1);
topPairs = rankedPairs(1:topN,:);
sigma = 1.2; % Smoothing factor for RF maps

% Preallocate
rf_overlap_on   = zeros(nPairsToRun,1);
rf_overlap_off  = zeros(nPairsToRun,1);
centroid_dist_on   = zeros(nPairsToRun,1); 
centroid_dist_off   = zeros(nPairsToRun,1);

for k = 1:nPairsToRun
    uA = rankedPairs(k,1);
    uB = rankedPairs(k,2);
    
    % 1. Process RFs: Mean across time, subtract baseline, and Smooth
    % Smoothing helps the thresholding identify 'blobs' rather than single pixels
    rfA_on  = imgaussfilt(mean(RFmap{uA}.ON.OnSet,3)  - RFmap{uA}.baseline, sigma);
    rfB_on  = imgaussfilt(mean(RFmap{uB}.ON.OnSet,3)  - RFmap{uB}.baseline, sigma);
    rfA_off = imgaussfilt(mean(RFmap{uA}.OFF.OnSet,3) - RFmap{uA}.baseline, sigma);
    rfB_off = imgaussfilt(mean(RFmap{uB}.OFF.OnSet,3) - RFmap{uB}.baseline, sigma);
    
    % 2. Thresholding (Top 20% of signal)
    threshA_on = prctile(rfA_on(:), 80);
    threshB_on = prctile(rfB_on(:), 80);
    threshA_off = prctile(rfA_off(:), 80);
    threshB_off = prctile(rfB_off(:), 80);
    
    maskA_on = rfA_on >= threshA_on;
    maskB_on = rfB_on >= threshB_on;
    maskA_off = rfA_off >= threshA_off;
    maskB_off = rfB_off >= threshB_off;

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
        rf_overlap_off(k) = intersection / union_area_off;
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

%% make scatter for all units (rf overlap vs coupling strength)

validIdx = ~isnan(rf_overlap_on) & ~isnan(coupling_wake_sorted);

x_data = rf_overlap_on(validIdx);
x_data_off = rf_overlap_off(validIdx);
y_data = coupling_wake_sorted(validIdx);
y_data_off = coupling_wake_sorted(validIdx);

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
[r, p] = corr(x_data_off, y_data_off);

figure('Color', 'w', 'Position', [100, 100, 800, 800]);
% Use 'binscatter' if you have many points, or 'scatter' with transparency
scatter(x_data_off, y_data_off, 20, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0.85 0.33 0.1]);
hold on;

% Add a trend line (Linear Regression)
coeffs = polyfit(x_data_off, y_data_off, 1);
fitX = linspace(min(x_data_off), max(x_data_off), 100);
fitY = polyval(coeffs, fitX);
plot(fitX, fitY, 'k', 'LineWidth', 1);

xlabel('OFF RF Overlap (Jaccard Index)');
ylabel('Wake Coupling Strength');
title(sprintf('Functional Coupling vs. Spatial Overlap\nr = %.3f (p = %.2e)', r, p));
grid on;


% Example for the Top 1 Pair
k = 1; 
uA = topPairs(k,1); uB = topPairs(k,2);

% Re-generate masks for visualization
m1 = maskA_on;
m2 = maskB_on;

% Create RGB image [Height x Width x 3]
sz = size(m1);
overlayImg = zeros(sz(1), sz(2), 3);
overlayImg(:,:,1) = m1; % Red channel
overlayImg(:,:,2) = m2; % Green channel

figure('Color', 'w');
imshow(overlayImg); 
title(sprintf('RF Overlap Pair %d (Cluster %d & %d)\nYellow = Shared Space', ...
    k, good_clusters(uA), good_clusters(uB)));
axis on;








% topPairs = rankedPairs(1:topN,:);
% 
% rf_overlap_on  = zeros(topN,1);
% rf_overlap_off = zeros(topN,1);
% 
% for k = 1:topN
% 
%     uA = topPairs(k,1);
%     uB = topPairs(k,2);
% 
%     % Collapse RF over time using mean
%     rfA_on  = mean(RFmap{uA}.ON.OnSet,3)  - RFmap{uA}.baseline;
%     rfB_on  = mean(RFmap{uB}.ON.OnSet,3)  - RFmap{uB}.baseline;
% 
%     rfA_off = mean(RFmap{uA}.OFF.OnSet,3) - RFmap{uA}.baseline;
%     rfB_off = mean(RFmap{uB}.OFF.OnSet,3) - RFmap{uB}.baseline;
% 
%     % 80th percentile threshold (keep top 20%)
%     threshA_on  = prctile(abs(rfA_on(:)), 80);
%     threshB_on  = prctile(abs(rfB_on(:)), 80);
% 
%     threshA_off = prctile(abs(rfA_off(:)), 80);
%     threshB_off = prctile(abs(rfB_off(:)), 80);
% 
%     % Binary masks
%     maskA_on  = abs(rfA_on)  >= threshA_on;
%     maskB_on  = abs(rfB_on)  >= threshB_on;
% 
%     maskA_off = abs(rfA_off) >= threshA_off;
%     maskB_off = abs(rfB_off) >= threshB_off;
% 
%     % Jaccard overlap
%     rf_overlap_on(k) = sum(maskA_on(:) & maskB_on(:)) / ...
%                        sum(maskA_on(:) | maskB_on(:));
% 
%     rf_overlap_off(k) = sum(maskA_off(:) & maskB_off(:)) / ...
%                         sum(maskA_off(:) | maskB_off(:));
% 
% end
% 
% figure;
% scatter(rf_overlap_on, coupling_wake_sorted(1:topN))
% xlabel('ON RF overlap (top 20%)')
% ylabel('Coupling strength (wake)')


























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





