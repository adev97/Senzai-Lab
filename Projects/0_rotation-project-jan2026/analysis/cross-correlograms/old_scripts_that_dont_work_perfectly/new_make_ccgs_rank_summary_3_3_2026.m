%% I am trying another time

%% -------------------------- SC CCG Analysis Script --------------------------
% Computes cross-correlograms (CCGs) for SC units, normalizes them using
% median baseline (excluding ±50 ms around 0), filters units by firing rate,
% and ranks synaptically coupled pairs for each state.

%% -------------------------- Setup Paths --------------------------
addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"));
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master"));
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\analysis"));

%% -------------------------- Load Kilosort Data --------------------------
sr = 30000;
ksDir = 'D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';

spikeTimes    = readNPY(fullfile(ksDir,'spike_times.npy')) + 1;
spikeClusters = readNPY(fullfile(ksDir,'spike_clusters.npy'));
templates     = readNPY(fullfile(ksDir,'templates.npy'));
chanPos       = readNPY(fullfile(ksDir,'channel_positions.npy')); % [chan x 2]

cluster_groups = readtable(fullfile(ksDir,'cluster_KSLabel.tsv'),'FileType','text','Delimiter','\t');
keepClusters   = cluster_groups.cluster_id(ismember(cluster_groups.KSLabel, {'good'}));

keepSpike = ismember(spikeClusters, keepClusters);
spikeTimes    = spikeTimes(keepSpike);
spikeClusters = spikeClusters(keepSpike);

%% -------------------------- Compute Unit Depth --------------------------
good_clusters = unique(spikeClusters);
nClusters = length(good_clusters);

unitDepth = zeros(nClusters,1);
for i = 1:nClusters
    clu = good_clusters(i);
    tempIdx = mode(spikeClusters(spikeClusters==clu));
    template = squeeze(templates(tempIdx+1,:,:));
    [~,peakChan] = max(max(abs(template),[],2));
    unitDepth(i) = chanPos(peakChan,2);
end

%% -------------------------- Map Spikes to Unit IDs --------------------------
ts = double(spikeTimes)/sr;
gs = double(spikeClusters);

unitIDs = zeros(size(gs));
for i = 1:nClusters
    unitIDs(gs == good_clusters(i)) = i;
end
s = [ts unitIDs];

%% -------------------------- Load Sleep States --------------------------
eegDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_eeg';
load(fullfile(eegDir, 'Mouse08_eeg.SleepState.states.mat'));

wakeInts = SleepState.ints.WAKEstate;
nremInts = SleepState.ints.NREMstate;
remInts  = SleepState.ints.REMstate;

s_wake = Restrict(s, wakeInts);
s_nrem = Restrict(s, nremInts);
s_rem  = Restrict(s, remInts);

%% -------------------------- Compute Firing Rates --------------------------
stateDurations = [sum(wakeInts(:,2)-wakeInts(:,1)), ...
                  sum(nremInts(:,2)-nremInts(:,1)), ...
                  sum(remInts(:,2)-remInts(:,1))];
              
all_s = {s_wake, s_nrem, s_rem};
minFiringRate = 0.5; % Hz

firingRates = zeros(nClusters,3);
for k = 1:nClusters
    firingRates(k,1) = sum(s_wake(:,2) == k)/stateDurations(1);
    firingRates(k,2) = sum(s_nrem(:,2) == k)/stateDurations(2);
    firingRates(k,3) = sum(s_rem(:,2) == k)/stateDurations(3);
end

fprintf('Units below %.2f Hz:\n', minFiringRate);
fprintf('  WAKE: %d / %d\n', sum(firingRates(:,1)<minFiringRate), nClusters);
fprintf('  NREM: %d / %d\n', sum(firingRates(:,2)<minFiringRate), nClusters);
fprintf('  REM : %d / %d\n', sum(firingRates(:,3)<minFiringRate), nClusters);

%% -------------------------- Compute and Normalize CCGs --------------------------
binSize      = 0.001; % 1 ms
duration     = 4;     % 8 s window
excludeWindow = 0.05; % ±50 ms for baseline
peakWindow    = 0.005; % ±5 ms for synaptic peak
states        = {'WAKE','NREM','REM'};
nStates       = length(states);

normCCGs       = cell(nStates,1);
peakVals       = cell(nStates,1);
validUnits     = cell(nStates,1);
validClusterIDs = cell(nStates,1);
rankedPairs    = cell(nStates,1);

for st = 1:nStates
    s_cur = all_s{st};
    stateDuration = stateDurations(st);
    
    % ------------------ Filter units by firing rate ------------------
    fr = zeros(nClusters,1);
    for k = 1:nClusters
        fr(k) = sum(s_cur(:,2) == k)/stateDuration;
    end
    vu = find(fr >= minFiringRate);                 % valid unit indices
    validUnits{st} = vu;
    validClusterIDs{st} = good_clusters(vu);       % original cluster IDs
    
    if length(vu)<2
        continue
    end
    
    % ------------------ Compute CCG ------------------
    [ccg, t] = CCG(s_cur(:,1), s_cur(:,2), 'binSize', binSize, 'duration', duration);
    
    baselineIdx = abs(t) > excludeWindow;
    peakIdx     = abs(t) <= peakWindow;
    
    nUnits = length(vu);
    normCCG = zeros(length(t), nUnits, nUnits);
    peakVal = zeros(nUnits, nUnits);
    
    % Only compute upper triangle to avoid duplicate pairs
    for i = 1:nUnits
        for j = i+1:nUnits
            raw = squeeze(ccg(:, vu(i), vu(j)));
            baseline = median(raw(baselineIdx));
            if baseline <= 0
                baseline = 1;
            end
            normCCG(:,i,j) = raw ./ baseline;
            peakVal(i,j) = max(normCCG(peakIdx,i,j));
        end
    end
    
    normCCGs{st} = normCCG;
    peakVals{st} = peakVal;
    
    fprintf('%s: %d valid units, %d unique pairs\n', states{st}, nUnits, nUnits*(nUnits-1)/2);
    
    % ------------------ Rank pairs ------------------
    pv = triu(peakVal,1);
    [vals, idx] = sort(pv(:),'descend');
    
    preID  = [];
    postID = [];
    peakStrength = [];
    
    for k = 1:length(idx)
        if vals(k) <= 1
            break
        end
        [i,j] = ind2sub(size(pv), idx(k));
        preID(end+1,1)  = validClusterIDs{st}(i);
        postID(end+1,1) = validClusterIDs{st}(j);
        peakStrength(end+1,1) = vals(k);
    end
    
    T = table(preID, postID, peakStrength);
    rankedPairs{st} = T;
    
    fprintf('%s: %d suprathreshold pairs found\n', states{st}, height(T));
    disp(T(1:min(10,height(T)),:))
end

%% -------------------------- Example: Plot top NREM pair --------------------------
% if ~isempty(rankedPairs{2}) % NREM = state 2
%     topPair = rankedPairs{2}(1,:);
%     i_idx = find(validClusterIDs{2} == topPair.preID);
%     j_idx = find(validClusterIDs{2} == topPair.postID);
% 
%     figure;
%     plot(t, squeeze(normCCGs{2}(:,i_idx,j_idx)));
%     xlabel('Time (s)');
%     ylabel('Normalized CCG');
%     title(sprintf('Top NREM pair: Cluster %d → Cluster %d', topPair.preID, topPair.postID));
% end

%% -------------------------- Full-matrix scatter: WAKE vs NREM --------------------------

%% Get normalized CCGs for WAKE and NREM
wakeCCG = normCCGs{1};  % [time x nWakeUnits x nWakeUnits]
nremCCG = normCCGs{2};  % [time x nNremUnits x nNremUnits]

% Valid units per state (indices into good_clusters)
unitsWake = validUnits{1};
unitsNREM = validUnits{2};

% Find units common to both states (by cluster ID)
commonClusterIDs = intersect(good_clusters(unitsWake), good_clusters(unitsNREM));

if isempty(commonClusterIDs)
    warning('No common units between WAKE and NREM.');
else
    % Map cluster IDs to indices in the normCCG arrays
    [~, wakeIdx] = ismember(commonClusterIDs, good_clusters(unitsWake));
    [~, nremIdx] = ismember(commonClusterIDs, good_clusters(unitsNREM));

    nUnits = length(commonClusterIDs);
    peakWindow = 0.005; % ±5 ms around 0
    peakIdx = abs(t) <= peakWindow;

    wakePeaks = [];
    nremPeaks = [];

    % Loop over all unique pairs (i < j) to avoid duplicates
    for i = 1:nUnits-1
        for j = i+1:nUnits
            wPeak = max(squeeze(wakeCCG(peakIdx, wakeIdx(i), wakeIdx(j))));
            nPeak = max(squeeze(nremCCG(peakIdx, nremIdx(i), nremIdx(j))));

            wakePeaks = [wakePeaks; wPeak];
            nremPeaks = [nremPeaks; nPeak];
        end
    end

    % Scatter plot
    figure;
    scatter(wakePeaks, nremPeaks, 50, 'filled');
    hold on;
    maxVal = max([wakePeaks; nremPeaks]);
    plot([0 maxVal], [0 maxVal], 'k--'); % unity line
    xlabel('WAKE coupling (norm CCG peak)');
    ylabel('NREM coupling (norm CCG peak)');
    title('SC unit pair coupling: WAKE vs NREM (all pairs)');
    axis square;
    grid on;

    % Correlation coefficient
    [R,P] = corrcoef(wakePeaks, nremPeaks);
    text(0.1*maxVal, 0.9*maxVal, sprintf('R = %.2f, p = %.3f', R(1,2), P(1,2)), 'FontSize', 12);

    % Show number of units/pairs
    fprintf('Scatter includes %d units forming %d unique pairs.\n', nUnits, length(wakePeaks));
end

%% -------------------------- Full-matrix scatter: WAKE vs REM --------------------------

% Get normalized CCGs for WAKE and REM
wakeCCG = normCCGs{1};  % [time x nWakeUnits x nWakeUnits]
remCCG  = normCCGs{3};  % [time x nREMUnits x nREMUnits]

% Valid units per state (indices into good_clusters)
unitsWake = validUnits{1};
unitsREM  = validUnits{3};

% Find units common to both states (by cluster ID)
commonClusterIDs = intersect(good_clusters(unitsWake), good_clusters(unitsREM));

if isempty(commonClusterIDs)
    warning('No common units between WAKE and REM.');
else
    % Map cluster IDs to indices in the normCCG arrays
    [~, wakeIdx] = ismember(commonClusterIDs, good_clusters(unitsWake));
    [~, remIdx]  = ismember(commonClusterIDs, good_clusters(unitsREM));

    nUnits = length(commonClusterIDs);
    peakWindow = 0.005; % ±5 ms around 0
    peakIdx = abs(t) <= peakWindow;

    wakePeaks = [];
    remPeaks  = [];

    % Loop over all unique pairs (i < j) to avoid duplicates
    for i = 1:nUnits-1
        for j = i+1:nUnits
            wPeak = max(squeeze(wakeCCG(peakIdx, wakeIdx(i), wakeIdx(j))));
            rPeak = max(squeeze(remCCG(peakIdx, remIdx(i), remIdx(j))));

            wakePeaks = [wakePeaks; wPeak];
            remPeaks  = [remPeaks; rPeak];
        end
    end

    % Scatter plot
    figure;
    scatter(wakePeaks, remPeaks, 50, 'filled');
    hold on;
    maxVal = max([wakePeaks; remPeaks]);
    plot([0 maxVal], [0 maxVal], 'k--'); % unity line
    xlabel('WAKE coupling (norm CCG peak)');
    ylabel('REM coupling (norm CCG peak)');
    title('SC unit pair coupling: WAKE vs REM (all pairs)');
    axis square;
    grid on;

    % Correlation coefficient
    [R,P] = corrcoef(wakePeaks, remPeaks);
    text(0.1*maxVal, 0.9*maxVal, sprintf('R = %.2f, p = %.3f', R(1,2), P(1,2)), 'FontSize', 12);

    % Show number of units/pairs
    fprintf('Scatter includes %d units forming %d unique pairs.\n', nUnits, length(wakePeaks));
end


%% Load external data for summary plotting (meanWav, RGC-SC-classification, RFMap)
meanWaveformDir   = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters';
classificationDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SC-RGC-classification-with-kmeans-all-channels-BEST';

load(fullfile(meanWaveformDir,   'meanWav_units.mat'));
load(fullfile(classificationDir, 'RGC-SC-classification.mat'));

load("\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_RFMapping\RF_maps\Spike_Data\Mouse08_20251007_810to2250_RFmap_SpikeRate.mat")

cell_type  = cell_type';
cluster_id = cluster_id';

%% -------------------------- Compute Receptive Field Overlap (Corrected Indexing) --------------------------

% Define states
stateLabels   = {'WAKE','NREM','REM'};
normCCG_all   = normCCGs;      % {1} = WAKE, {2} = NREM, {3} = REM
validUnits_all = validUnits;   % valid units per state (indices into good_clusters)

% Find units present in all three states
commonUnits = intersect(intersect(validUnits_all{1}, validUnits_all{2}), validUnits_all{3});
nCommon = length(commonUnits);

if nCommon < 2
    warning('Less than 2 common units across states. Skipping RF overlap computation.');
else
    % Generate all unique pairs of common units
    pairs_idx = nchoosek(1:nCommon, 2);  % [i j] pairs
    nPairs = size(pairs_idx,1);

    % Map commonUnits to indices in normCCGs for each state
    wakeMap = find(ismember(validUnits_all{1}, commonUnits));
    nremMap = find(ismember(validUnits_all{2}, commonUnits));
    remMap  = find(ismember(validUnits_all{3}, commonUnits));

    peakWindow = 0.005;
    peakIdx = abs(t) <= peakWindow;

    % Preallocate peak arrays
    wakePeaks = nan(nPairs,1);
    nremPeaks = nan(nPairs,1);
    remPeaks  = nan(nPairs,1);

    % Compute peaks for all pairs
    for p = 1:nPairs
        i = pairs_idx(p,1);
        j = pairs_idx(p,2);

        wakePeaks(p) = max(squeeze(normCCG_all{1}(peakIdx, wakeMap(i), wakeMap(j))));
        nremPeaks(p) = max(squeeze(normCCG_all{2}(peakIdx, nremMap(i), nremMap(j))));
        remPeaks(p)  = max(squeeze(normCCG_all{3}(peakIdx, remMap(i), remMap(j))));
    end

    fprintf('Computed CCG peaks for %d common pairs.\n', nPairs);

    %% -------------------------- Compute RF Overlap (Jaccard) --------------------------
    jaccard_ON  = nan(nPairs,1);
    jaccard_OFF = nan(nPairs,1);

    sigma = 1.2;
    hSize = 2*ceil(2*sigma)+1;
    h = fspecial('gaussian', hSize, sigma);
    maskVal = 80;

    for p = 1:nPairs
        i = pairs_idx(p,1);
        j = pairs_idx(p,2);

        unitA = commonUnits(i);
        unitB = commonUnits(j);

        % Smooth RFs
        rfA_on  = imfilter(mean(RFmap{unitA}.ON.OnSet,3)  - RFmap{unitA}.baseline, h,'replicate');
        rfB_on  = imfilter(mean(RFmap{unitB}.ON.OnSet,3)  - RFmap{unitB}.baseline, h,'replicate');
        rfA_off = imfilter(mean(RFmap{unitA}.OFF.OnSet,3) - RFmap{unitA}.baseline, h,'replicate');
        rfB_off = imfilter(mean(RFmap{unitB}.OFF.OnSet,3) - RFmap{unitB}.baseline, h,'replicate');

        noise_A_on  = 2*std(rfA_on(:));
        noise_A_off = 2*std(rfA_off(:));
        noise_B_on  = 2*std(rfB_on(:));
        noise_B_off = 2*std(rfB_off(:));

        % ON masks
        mA_on = (rfA_on >= prctile(rfA_on(:), maskVal)) & rfA_on > noise_A_on;
        mB_on = (rfB_on >= prctile(rfB_on(:), maskVal)) & rfB_on > noise_B_on;

        % OFF masks
        mA_off = (rfA_off >= prctile(rfA_off(:), maskVal)) & rfA_off > noise_A_off;
        mB_off = (rfB_off >= prctile(rfB_off(:), maskVal)) & rfB_off > noise_B_off;

        % Compute Jaccard indices
        jaccard_ON(p)  = computeRFJaccard(mA_on,  mB_on);
        jaccard_OFF(p) = computeRFJaccard(mA_off, mB_off);
    end

    fprintf('Computed RF overlap (Jaccard) for all %d pairs.\n', nPairs);

    %% -------------------------- Plot RF Overlap vs Coupling --------------------------
    figure('Color','w','Position',[100,100,1400,900]);
    couplingAll = {wakePeaks, nremPeaks, remPeaks};
    stateColors = {[0 0 0], [0.8 0.2 0.2], [0.2 0.4 0.8]};

    for col = 1:3
        coupling = couplingAll{col};

        % --- ON row ---
        subplot(2,3,col);
        validIdx = ~isnan(jaccard_ON) & ~isnan(coupling);
        x = jaccard_ON(validIdx);
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

        % --- OFF row ---
        subplot(2,3,col+3);
        validIdx = ~isnan(jaccard_OFF) & ~isnan(coupling);
        x = jaccard_OFF(validIdx);
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
end

% %% -------------------------- Compute Receptive Field Overlap --------------------------
% 
% %% 1. Jaccard Index
% wakeCCG = normCCGs{1};   % make explicit
% nremCCG = normCCGs{2};
% remCCG  = normCCGs{3};
% 
% unitsWake = validUnits{1};              % indices into good_clusters
% unitsNREM = validUnits{2};
% unitsREM  = validUnits{3};
% 
% commonUnits = intersect(intersect(unitsWake, unitsNREM), unitsREM);
% 
% wakeClusterIDs = good_clusters(unitsWake);
% 
% % generate pairs
% nUnits = length(commonUnits);
% pairs_idx = zeros(nchoosek(nUnits,2),2);
% k=1;
% for i = 1:nUnits-1
%     for j = i+1:nUnits
%         pairs_idx(k,:) = [i j];
%         k = k + 1;
%     end
% end
% 
% nPairs = size(pairs_idx,1);
% 
% peakWindow = 0.005;
% peakIdx = abs(t) <= peakWindow;
% 
% wakePeaks = nan(nPairs,1);
% nremPeaks = nan(nPairs,1);
% remPeaks  = nan(nPairs,1);
% 
% for p = 1:nPairs
%     i = pairs_idx(p,1);
%     j = pairs_idx(p,2);
% 
%     wakePeaks(p) = max(squeeze(normCCGs{1}(peakIdx, i, j)));
%     nremPeaks(p) = max(squeeze(normCCGs{2}(peakIdx, i, j)));
%     remPeaks(p)  = max(squeeze(normCCGs{3}(peakIdx, i, j)));
% end
% 
% 
% 
% % wakePeaks = nan(nPairs,1);
% % nremPeaks = nan(nPairs,1);
% % remPeaks  = nan(nPairs,1);
% % 
% % peakWindow = 0.005;        % ±5 ms
% % peakIdx = abs(t) <= peakWindow;
% % 
% % for p = 1:nPairs
% %     i = pairs_idx(p,1);
% %     j = pairs_idx(p,2);
% % 
% %     wakePeaks(p) = max(squeeze(normCCGs{1}(peakIdx, i, j)));
% %     nremPeaks(p) = max(squeeze(normCCGs{2}(peakIdx, i, j)));
% %     remPeaks(p)  = max(squeeze(normCCGs{3}(peakIdx, i, j)));
% % end
% 
% fprintf('Total common pairs: %d\n', nPairs);
% 
% % extract wake coupling for all pairs
% % peakWindow = 0.005;
% % peakIdx = abs(t) <= peakWindow;
% % 
% % wakeCoupling = nan(nPairs,1);
% % 
% % for p = 1:pairCount
% %     i = pairs_idx(p,1);
% %     j = pairs_idx(p,2);
% % 
% %     wakeCoupling(p) = max(squeeze(wakeCCG(peakIdx,i,j)));
% % end
% 
% % compute jaccard
% jaccard_ON  = nan(nPairs,1);
% jaccard_OFF = nan(nPairs,1);
% 
% sigma = 1.2;
% hSize = 2*ceil(2*sigma)+1;
% h = fspecial('gaussian',hSize,sigma);
% 
% maskVal = 80;
% 
% for p = 1:nPairs
% 
%     i = pairs_idx(p,1);
%     j = pairs_idx(p,2);
% 
%     unitA = commonUnits(i);   % index into RFmap
%     unitB = commonUnits(j);
% 
%     % ---- Smooth RFs ----
%     rfA_on  = imfilter(mean(RFmap{unitA}.ON.OnSet,3)  - RFmap{unitA}.baseline,h,'replicate');
%     rfB_on  = imfilter(mean(RFmap{unitB}.ON.OnSet,3)  - RFmap{unitB}.baseline,h,'replicate');
%     rfA_off = imfilter(mean(RFmap{unitA}.OFF.OnSet,3) - RFmap{unitA}.baseline,h,'replicate');
%     rfB_off = imfilter(mean(RFmap{unitB}.OFF.OnSet,3) - RFmap{unitB}.baseline,h,'replicate');
% 
%     noise_A_on  = 2*std(rfA_on(:));
%     noise_A_off = 2*std(rfA_off(:));
%     noise_B_on  = 2*std(rfB_on(:));
%     noise_B_off = 2*std(rfB_off(:));
% 
%     % ON masks
%     mA_on = (rfA_on >= prctile(rfA_on(:),maskVal)) & rfA_on > noise_A_on;
%     mB_on = (rfB_on >= prctile(rfB_on(:),maskVal)) & rfB_on > noise_B_on;
% 
%     % OFF masks
%     mA_off = (rfA_off >= prctile(rfA_off(:),maskVal)) & rfA_off > noise_A_off;
%     mB_off = (rfB_off >= prctile(rfB_off(:),maskVal)) & rfB_off > noise_B_off;
% 
%     % Compute separate Jaccards
%     jaccard_ON(p)  = computeRFJaccard(mA_on,  mB_on);
%     jaccard_OFF(p) = computeRFJaccard(mA_off, mB_off);
% 
%     % Optional combined
%     % jaccard_combined(p) = computeRFJaccard(mA_on | mA_off, mB_on | mB_off);
% end
% 
% fprintf('Computed RF overlap for all pairs.\n');
% 
% % plot rf overlap vs coupling strength wake
% % figure;
% % scatter(jaccard_combined, wakeCoupling, 30, 'filled');
% % xlabel('RF overlap (Jaccard)');
% % ylabel('WAKE coupling');
% % axis square;
% % grid on;
% % 
% % [R,P] = corrcoef(jaccard_combined, wakeCoupling,'Rows','complete');
% % title(sprintf('Wake Coupling vs Jaccard Overlap Index\n R = %.2f, p = %.3g',R(1,2),P(1,2)));
% % lsline;   % least squares regression line
% 
% 
% 
% %% Try with log scale
% % semilogx(jaccard_combined, wakeCoupling,'.')
% 
% %% -------------------------- Plot RF Overlap vs Coupling --------------------------
% 
% figure('Color', 'w', 'Position', [100, 100, 1400, 900]);
% 
% stateLabels   = {'WAKE', 'NREM', 'REM'};
% couplingAll   = {wakePeaks, nremPeaks, remPeaks};
% stateColors   = {[0 0 0], [0.8 0.2 0.2], [0.2 0.4 0.8]};
% 
% rf_overlap_on  = jaccard_ON;
% rf_overlap_off = jaccard_OFF;
% 
% for col = 1:3  % columns = states (WAKE, NREM, REM)
% 
%     coupling = couplingAll{col};
% 
%     % --- ON row (row 1) ---
%     subplot(2,3,col);
%     validIdx = ~isnan(rf_overlap_on) & ~isnan(coupling);
%     x = rf_overlap_on(validIdx);
%     y = coupling(validIdx);
% 
%     [r, p] = corr(x, y);
%     scatter(x, y, 15, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', stateColors{col});
%     hold on;
%     coeffs = polyfit(x, y, 1);
%     fitX = linspace(min(x), max(x), 100);
%     plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1.2);
%     xlabel('ON RF Overlap (Jaccard)');
%     ylabel([stateLabels{col} ' Coupling']);
%     title(sprintf('%s — ON\nr = %.3f, p = %.2e', stateLabels{col}, r, p));
%     grid on;
% 
%     % --- OFF row (row 2) ---
%     subplot(2,3,col+3);
%     validIdx = ~isnan(rf_overlap_off) & ~isnan(coupling);
%     x = rf_overlap_off(validIdx);
%     y = coupling(validIdx);
% 
%     [r, p] = corr(x, y);
%     scatter(x, y, 15, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', stateColors{col});
%     hold on;
%     coeffs = polyfit(x, y, 1);
%     fitX = linspace(min(x), max(x), 100);
%     plot(fitX, polyval(coeffs, fitX), 'k', 'LineWidth', 1.2);
%     xlabel('OFF RF Overlap (Jaccard)');
%     ylabel([stateLabels{col} ' Coupling']);
%     title(sprintf('%s — OFF\nr = %.3f, p = %.2e', stateLabels{col}, r, p));
%     grid on;
% 
% end
% 
% sgtitle('RF Overlap vs Coupling Strength Across States', 'FontSize', 16, 'FontWeight', 'bold');

% nPairs = size(pairs,1);
% 
% jaccard_ON  = nan(nPairs,1);
% jaccard_OFF = nan(nPairs,1);
% jaccard_combined = nan(nPairs,1);
% 
% sigma = 1.2;
% hSize = 2*ceil(2*sigma)+1;
% h = fspecial('gaussian',hSize,sigma);
% 
% for p = 1:nPairs
% 
%     unitA = pairs(p,1);
%     unitB = pairs(p,2);
% 
%     % ---- Smooth RFs ----
%     rfA_on  = imfilter(mean(RFmap{unitA}.ON.OnSet,3)  - RFmap{unitA}.baseline,h,'replicate');
%     rfB_on  = imfilter(mean(RFmap{unitB}.ON.OnSet,3)  - RFmap{unitB}.baseline,h,'replicate');
%     rfA_off = imfilter(mean(RFmap{unitA}.OFF.OnSet,3) - RFmap{unitA}.baseline,h,'replicate');
%     rfB_off = imfilter(mean(RFmap{unitB}.OFF.OnSet,3) - RFmap{unitB}.baseline,h,'replicate');
% 
%     % ---- Threshold masks ----
%     noise_A_on  = 2*std(rfA_on(:));
%     noise_B_on  = 2*std(rfB_on(:));
%     noise_A_off = 2*std(rfA_off(:));
%     noise_B_off = 2*std(rfB_off(:));
% 
%     mA_on  = (rfA_on  >= prctile(rfA_on(:),maskVal))  & (rfA_on  > noise_A_on);
%     mB_on  = (rfB_on  >= prctile(rfB_on(:),maskVal))  & (rfB_on  > noise_B_on);
%     mA_off = (rfA_off >= prctile(rfA_off(:),maskVal)) & (rfA_off > noise_A_off);
%     mB_off = (rfB_off >= prctile(rfB_off(:),maskVal)) & (rfB_off > noise_B_off);
% 
%     % ---- Compute Jaccard ----
%     jaccard_ON(p)  = computeRFJaccard(mA_on,  mB_on);
%     jaccard_OFF(p) = computeRFJaccard(mA_off, mB_off);
% 
%     % Combined ON+OFF mask
%     mA_comb = mA_on | mA_off;
%     mB_comb = mB_on | mB_off;
%     jaccard_combined(p) = computeRFJaccard(mA_comb, mB_comb);
% 
% end
% 
% fprintf('Computed RF overlap for %d pairs.\n', nPairs);


%% ------------------ Make summary figures for highest ranked pairs units ---------------------

pngOutputDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\00000000-PAIR-SUMMARIES-WITH-RF-3-4-2026';
if ~exist(pngOutputDir,'dir'), mkdir(pngOutputDir); end

all_ccgs = {normCCGs{1}, normCCGs{2}, normCCGs{3}}; % WAKE, NREM, REM
% t        = t;                                        % time vector from CCG
pairs    = [rankedPairs{1}.preID, rankedPairs{1}.postID];  % WAKE top pairs (cluster IDs)

% Map cluster IDs to indices in good_clusters
[~, pairIdxA] = ismember(pairs(:,1), good_clusters);
[~, pairIdxB] = ismember(pairs(:,2), good_clusters);
pairs_idx = [pairIdxA, pairIdxB];  % indices into your arrays

% Coupling values for each state
coupling_wake = wakePeaks;   % you already computed these
coupling_nrem = nremPeaks;
coupling_rem  = remPeaks;

% Other inputs (from your workspace)
sortBy = 'wake';  % sort top pairs by WAKE coupling
topN   = 10;      % top 30 pairs
pdfDir = 'D:\SummaryFigures\Top30';  % folder to save PNGs
maskVal = 90;     % percentile threshold for RF masks
dispTime = 0.5;   % ±0.5 s around 0 for CCG plots

% Call the function
makePairSummaryWithRFsTHISONEWORKS(all_ccgs, t, pairs_idx, good_clusters, ...
                       coupling_wake, coupling_nrem, coupling_rem, ...
                       unitDepth, meanWav, cell_type, ...
                       xpos, ypos, sr, sortBy, topN, pngOutputDir, RFmap, maskVal, dispTime);




