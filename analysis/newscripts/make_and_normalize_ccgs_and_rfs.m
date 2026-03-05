% correctly_normalized_ccgs

%% Start from Scratch, 3/4/2026

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

%% ----------- Calculate CCGs -----------------------

binSize  = 0.001;   % 1 ms bins
duration = 0.1;     % ±50 ms window (good for spike sorting / monosynaptic)

spikeTimes_wake = cell(nClusters, 1);
spikeTimes_nrem = cell(nClusters, 1);
spikeTimes_rem  = cell(nClusters, 1);

for k = 1:nClusters
    spikeTimes_wake{k} = s_wake(s_wake(:,2) == k, 1);
    spikeTimes_nrem{k} = s_nrem(s_nrem(:,2) == k, 1);
    spikeTimes_rem{k}  = s_rem(s_rem(:,2)  == k, 1);
end

fprintf('Computing CCGs...\n');

[ccg_wake, t_ccg] = CCG(spikeTimes_wake, [], 'binSize', binSize, 'duration', duration);
[ccg_nrem, ~] = CCG(spikeTimes_nrem, [], 'binSize', binSize, 'duration', duration);
[ccg_rem,  ~] = CCG(spikeTimes_rem,  [], 'binSize', binSize, 'duration', duration);

fprintf('  CCGs computed.');

%% Save CCGs
saveFile = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\000000000000-final_figures\saved_ccgs\100ms';
if ~exist(saveFile, 'dir'); mkdir(saveFile); end

save(fullfile(saveFile, 'ccgs_all.mat'), 'ccg_wake', 'ccg_nrem', 'ccg_rem', 't_ccg');

fprintf(' Begin Normalization...\n');

%% ---------- Normalize CCGs and filter for low spike count -----------

%Precompute reference spike counts for normalization
nRefSpikes_wake = zeros(nClusters, 1);
nRefSpikes_nrem = zeros(nClusters, 1);
nRefSpikes_rem  = zeros(nClusters, 1);

for k = 1:nClusters
    nRefSpikes_wake(k) = sum(s_wake(:,2) == k);
    nRefSpikes_nrem(k) = sum(s_nrem(:,2) == k);
    nRefSpikes_rem(k)  = sum(s_rem(:,2)  == k);
end

% Normalize to rate
ccg_wake_rate = double(ccg_wake) ./ reshape(nRefSpikes_wake, [1 nClusters 1]) ./ binSize;
ccg_nrem_rate = double(ccg_nrem) ./ reshape(nRefSpikes_nrem, [1 nClusters 1]) ./ binSize;
ccg_rem_rate  = double(ccg_rem)  ./ reshape(nRefSpikes_rem,  [1 nClusters 1]) ./ binSize;

lowRate_wake = firingRate_wake < minFiringRate;
lowRate_nrem = firingRate_nrem < minFiringRate;
lowRate_rem  = firingRate_rem  < minFiringRate;

ccg_wake_rate(:, lowRate_wake, :) = NaN;  ccg_wake_rate(:, :, lowRate_wake) = NaN;
ccg_nrem_rate(:, lowRate_nrem, :) = NaN;  ccg_nrem_rate(:, :, lowRate_nrem) = NaN;
ccg_rem_rate(:,  lowRate_rem,  :) = NaN;  ccg_rem_rate(:,  :, lowRate_rem)  = NaN;

fprintf('Valid units after firing rate threshold (%.1f Hz):\n', minFiringRate);
fprintf('  WAKE: %d / %d\n', sum(~lowRate_wake), nClusters);
fprintf('  NREM: %d / %d\n', sum(~lowRate_nrem), nClusters);
fprintf('  REM : %d / %d\n', sum(~lowRate_rem),  nClusters);


%% ------------- plot example unit ccg (raw rates, not counts) ------------
ksIDi = 370; % kilosort cluster id
ksIDj = 376; % kilosort cluster id

ui = find(good_clusters == ksIDi);  % reference unit
uj = find(good_clusters == ksIDj);  % target unit

figure;
hold on;
plot(t_ccg * 1000, ccg_wake_rate(:, ui, uj), 'k', 'LineWidth', 1.5);
plot(t_ccg * 1000, ccg_nrem_rate(:, ui, uj), 'r', 'LineWidth', 1.5);
plot(t_ccg * 1000, ccg_rem_rate(:,  ui, uj), 'b', 'LineWidth', 1.5);

xlabel('Time lag (ms)');
ylabel('Firing rate (spk/s)');
title(sprintf('CCG: %d → %d  (units %d, %d)', ksIDi, ksIDj, ui, uj));
legend('WAKE', 'NREM', 'REM', 'Location', 'best');
xlim([t_ccg(1)*1000, t_ccg(end)*1000]);
box off;
hold off;


%% ------------ Rank pairs
% monoWin  = 0.005;
% monoIdx  = abs(t_ccg) <= monoWin;
% flankIdx = abs(t_ccg) >= 0.8 * (duration/2);
% 
% monoScore_wake = zeros(nClusters, nClusters);
% monoScore_nrem = zeros(nClusters, nClusters);
% monoScore_rem  = zeros(nClusters, nClusters);
% 
% for ui = 1:nClusters
%     for uj = 1:nClusters
%         if ui == uj; continue; end
% 
%         wake_trace = ccg_wake_rate(:, ui, uj);
%         nrem_trace = ccg_nrem_rate(:, ui, uj);
%         rem_trace  = ccg_rem_rate(:,  ui, uj);
% 
%         if ~all(isnan(wake_trace))
%             monoScore_wake(ui,uj) = max(wake_trace(monoIdx)) - mean(wake_trace(flankIdx), 'omitnan');
%         end
%         if ~all(isnan(nrem_trace))
%             monoScore_nrem(ui,uj) = max(nrem_trace(monoIdx)) - mean(nrem_trace(flankIdx), 'omitnan');
%         end
%         if ~all(isnan(rem_trace))
%             monoScore_rem(ui,uj)  = max(rem_trace(monoIdx))  - mean(rem_trace(flankIdx), 'omitnan');
%         end
%     end
% end
% 
% [ii, jj]  = ndgrid(1:nClusters, 1:nClusters);
% mask      = ii ~= jj;
% ui_list   = ii(mask);
% uj_list   = jj(mask);
% 
% ranked_wake = sortrows(table(ui_list, uj_list, good_clusters(ui_list), good_clusters(uj_list), monoScore_wake(mask), ...
%     'VariableNames', {'ui','uj','ksID_ref','ksID_target','monoScore'}), 'monoScore', 'descend');
% ranked_nrem = sortrows(table(ui_list, uj_list, good_clusters(ui_list), good_clusters(uj_list), monoScore_nrem(mask), ...
%     'VariableNames', {'ui','uj','ksID_ref','ksID_target','monoScore'}), 'monoScore', 'descend');
% ranked_rem  = sortrows(table(ui_list, uj_list, good_clusters(ui_list), good_clusters(uj_list), monoScore_rem(mask),  ...
%     'VariableNames', {'ui','uj','ksID_ref','ksID_target','monoScore'}), 'monoScore', 'descend');

%% ---------------------------- Normalize ALL CCGs pairs
monoWin  = 0.005;
monoIdx  = abs(t_ccg) <= monoWin;
flankIdx = abs(t_ccg) >= 0.8 * (duration/2);

minBaseline  = 1.0;
minRefSpikes = 50;

nBins = length(t_ccg);
ccg_wake_norm = nan(nBins, nClusters, nClusters);
ccg_nrem_norm = nan(nBins, nClusters, nClusters);
ccg_rem_norm  = nan(nBins, nClusters, nClusters);

for ui = 1:nClusters
    for uj = 1:nClusters
        if ui == uj; continue; end

        [ccg_wake_norm(:,ui,uj), ~, ~] = normalize_ccg(ccg_wake_rate(:,ui,uj), flankIdx, minBaseline, minRefSpikes, nRefSpikes_wake(ui));
        [ccg_nrem_norm(:,ui,uj), ~, ~] = normalize_ccg(ccg_nrem_rate(:,ui,uj), flankIdx, minBaseline, minRefSpikes, nRefSpikes_nrem(ui));
        [ccg_rem_norm(:,ui,uj),  ~, ~] = normalize_ccg(ccg_rem_rate(:, ui,uj), flankIdx, minBaseline, minRefSpikes, nRefSpikes_rem(ui));
    end
end

fprintf('Normalized CCG matrices built.\n');
fprintf('  Valid pairs - WAKE: %d  NREM: %d  REM: %d\n', ...
    sum(~isnan(squeeze(ccg_wake_norm(1,:,:))), 'all'), ...
    sum(~isnan(squeeze(ccg_nrem_norm(1,:,:))), 'all'), ...
    sum(~isnan(squeeze(ccg_rem_norm(1,:,:))),  'all'));

%% ---------------------- Rank pairs on normalized CCG
monoScore_wake = zeros(nClusters, nClusters);
monoScore_nrem = zeros(nClusters, nClusters);
monoScore_rem  = zeros(nClusters, nClusters);

for ui = 1:nClusters
    for uj = 1:nClusters
        if ui == uj; continue; end

        wake_trace = ccg_wake_norm(:, ui, uj);
        nrem_trace = ccg_nrem_norm(:, ui, uj);
        rem_trace  = ccg_rem_norm(:,  ui, uj);

        if ~all(isnan(wake_trace))
            monoScore_wake(ui,uj) = max(wake_trace(monoIdx)) - 1;
        end
        if ~all(isnan(nrem_trace))
            monoScore_nrem(ui,uj) = max(nrem_trace(monoIdx)) - 1;
        end
        if ~all(isnan(rem_trace))
            monoScore_rem(ui,uj)  = max(rem_trace(monoIdx))  - 1;
        end
    end
end

% Build ranked tables
[ii, jj]  = ndgrid(1:nClusters, 1:nClusters);
mask      = ii ~= jj;
ui_list   = ii(mask);
uj_list   = jj(mask);

ranked_wake = sortrows(table(ui_list, uj_list, good_clusters(ui_list), good_clusters(uj_list), monoScore_wake(mask), ...
    'VariableNames', {'ui','uj','ksID_ref','ksID_target','monoScore'}), 'monoScore', 'descend');
ranked_nrem = sortrows(table(ui_list, uj_list, good_clusters(ui_list), good_clusters(uj_list), monoScore_nrem(mask), ...
    'VariableNames', {'ui','uj','ksID_ref','ksID_target','monoScore'}), 'monoScore', 'descend');
ranked_rem  = sortrows(table(ui_list, uj_list, good_clusters(ui_list), good_clusters(uj_list), monoScore_rem(mask),  ...
    'VariableNames', {'ui','uj','ksID_ref','ksID_target','monoScore'}), 'monoScore', 'descend');

%% Remove duplicate pairs in ranking
% After building ranked_wake, remove duplicate pairs keeping higher score
seen = false(nClusters, nClusters);
keepRows = false(height(ranked_wake), 1);

for row = 1:height(ranked_wake)
    ui = ranked_wake.ui(row);
    uj = ranked_wake.uj(row);
    if ~seen(ui,uj) && ~seen(uj,ui)
        keepRows(row) = true;
        seen(ui,uj) = true;
        seen(uj,ui) = true;
    end
end

ranked_wake_dedup = ranked_wake(keepRows, :);

%% Plot top ranked pair valid across all states
for row = 1:height(ranked_wake)
    ui = ranked_wake.ui(row);
    uj = ranked_wake.uj(row);
    if ~all(isnan(ccg_nrem_norm(:,ui,uj))) && ~all(isnan(ccg_rem_norm(:,ui,uj)))
        break
    end
end

fprintf('Plotting top valid pair: KS%d -> KS%d (rank %d)\n', ...
    good_clusters(ui), good_clusters(uj), row);

wake_norm = ccg_wake_norm(:, ui, uj);
nrem_norm = ccg_nrem_norm(:, ui, uj);
rem_norm  = ccg_rem_norm(:,  ui, uj);

wake_ok = ~all(isnan(wake_norm));
nrem_ok = ~all(isnan(nrem_norm));
rem_ok  = ~all(isnan(rem_norm));

figure; hold on;
if wake_ok; plot(t_ccg*1000, wake_norm, 'k', 'LineWidth', 1.5); end
if nrem_ok; plot(t_ccg*1000, nrem_norm, 'r', 'LineWidth', 1.5); end
if rem_ok;  plot(t_ccg*1000, rem_norm,  'b', 'LineWidth', 1.5); end

yline(1, '--', 'Color', [0.5 0.5 0.5], 'Label', 'baseline');
xline(0, '--', 'Color', [0.5 0.5 0.5]);
xline([-5 5], ':', 'Color', [0.7 0.7 0.7]);
xlabel('Time lag (ms)');
ylabel('Normalized rate (a.u.)');
title(sprintf('CCG: KS%d → KS%d  (units %d, %d)', good_clusters(ui), good_clusters(uj), ui, uj));

legendEntries = {};
if wake_ok; legendEntries{end+1} = 'WAKE'; end
if nrem_ok; legendEntries{end+1} = 'NREM'; end
if rem_ok;  legendEntries{end+1} = 'REM';  end
legend(legendEntries, 'Location', 'best');
box off;

%% --------- Load External Items -----------

% load RFmaps
load("\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08_RFMapping\RF_maps\Spike_Data\Mouse08_20251007_810to2250_RFmap_SpikeRate.mat");

% Load meanWav and cell types
meanWaveformDir   = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters';
classificationDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SC-RGC-classification-with-kmeans-all-channels-BEST';

load(fullfile(meanWaveformDir,   'meanWav_units.mat'));
load(fullfile(classificationDir, 'RGC-SC-classification.mat'));

xpos = chanPos(:,1);
ypos = chanPos(:,2);

% Flip cell_type orientation if needed
cell_type  = cell_type';
cluster_id = cluster_id';


%% ---------- Make Summary Figures -------

% Load cell type classification
% Expected variables from .mat: cluster_id (array), cell_type (cell array of strings)
% Map cluster_id -> cell_type string for fast lookup
cellTypeMap = containers.Map(num2cell(double(cluster_id)), cell_type);

% Output directory
pngOutputDir = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\000000000000-final_figures";
saveDir = fullfile(pngOutputDir, 'dedup-100ms-CCG_SummaryFigures');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

nSummary = 30;
figCount = 0;

for row = 1:height(ranked_wake_dedup) % or ranked_wake
    ui = ranked_wake_dedup.ui(row);
    uj = ranked_wake_dedup.uj(row);

    % Check both units valid across all states
    if all(isnan(ccg_nrem_norm(:,ui,uj))) || all(isnan(ccg_rem_norm(:,ui,uj)))
        continue
    end

    figCount = figCount + 1;
    if figCount > nSummary; break; end

    ksIDi = good_clusters(ui);
    ksIDj = good_clusters(uj);

    % --- Get normalized CCG traces ---
    wake_tr = ccg_wake_norm(:, ui, uj);
    nrem_tr = ccg_nrem_norm(:, ui, uj);
    rem_tr  = ccg_rem_norm(:,  ui, uj);

    wake_ok = ~all(isnan(wake_tr));
    nrem_ok = ~all(isnan(nrem_tr));
    rem_ok  = ~all(isnan(rem_tr));

    % --- Get cell types ---
    if isKey(cellTypeMap, double(ksIDi))
        ctI = cellTypeMap(double(ksIDi));
    else
        ctI = 'unknown';
    end
    if isKey(cellTypeMap, double(ksIDj))
        ctJ = cellTypeMap(double(ksIDj));
    else
        ctJ = 'unknown';
    end

    % --- Get RF maps ---
    onI  = mean(RFmap{ui}.ON.OnSet,  3);  % [y x x]
    offI = mean(RFmap{ui}.OFF.OnSet, 3);
    onJ  = mean(RFmap{uj}.ON.OnSet,  3);
    offJ = mean(RFmap{uj}.OFF.OnSet, 3);

    % ---------------------------------------------------------------
    % Build figure
    % Layout:
    %   Row 1 (tall): CCG spanning full width
    %   Row 2: [waveform i] [waveform j]
    %   Row 3: [cell type i] [cell type j]
    %   Row 4: [ON RF i]    [ON RF j]
    %   Row 5: [OFF RF i]   [OFF RF j]
    % ---------------------------------------------------------------
    fig = figure('Units','centimeters','Position',[2 2 28 34],'Visible','on');

    % --- CCG (top, full width) ---
    ax_ccg = subplot(5, 2, [1 2]);
    hold(ax_ccg, 'on');
    if wake_ok; plot(ax_ccg, t_ccg*1000, wake_tr, 'k', 'LineWidth', 1.5); end
    if nrem_ok; plot(ax_ccg, t_ccg*1000, nrem_tr, 'r', 'LineWidth', 1.5); end
    if rem_ok;  plot(ax_ccg, t_ccg*1000, rem_tr,  'b', 'LineWidth', 1.5); end
    yline(ax_ccg, 1, '--', 'Color', [0.5 0.5 0.5]);
    xline(ax_ccg, 0, '--', 'Color', [0.5 0.5 0.5]);
    xline(ax_ccg, [-5 5], ':', 'Color', [0.7 0.7 0.7]);
    xlabel(ax_ccg, 'Time lag (ms)');
    ylabel(ax_ccg, 'Normalized rate');
    title(ax_ccg, sprintf('CCG: %d → %d  |  Rank %d  |  Score: %.2f', ...
        ksIDi, ksIDj, row, ranked_wake_dedup.monoScore(row)));
    legEntries = {};
    if wake_ok; legEntries{end+1} = 'WAKE'; end
    if nrem_ok; legEntries{end+1} = 'NREM'; end
    if rem_ok;  legEntries{end+1} = 'REM';  end
    legend(ax_ccg, legEntries, 'Location','best');
    box(ax_ccg, 'off');

    % --- Waveforms --- 
    
    % Unit i
    ax_wI = subplot(5, 2, 3);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, ui, ax_wI);
    
    % Unit j
    ax_wJ = subplot(5, 2, 4);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, uj, ax_wJ);

    % --- Cell type labels ---
    ax_ctI = subplot(5, 2, 5);
    axis(ax_ctI, 'off');
    text(0.5, 0.5, sprintf('Unit %d\n%s', ksIDi, ctI), 'Parent', ax_ctI, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', get_cell_type_color(ctI));


    ax_ctJ = subplot(5, 2, 6);
    axis(ax_ctJ, 'off');
    text(0.5, 0.5, sprintf('Unit %d\n%s', ksIDj, ctJ), 'Parent', ax_ctJ, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', get_cell_type_color(ctJ));

    % --- ON RF ---
    ax_onI = subplot(5, 2, 7);
    imagesc(ax_onI, onI);
    colormap(ax_onI, gray);
    axis(ax_onI, 'off');
    title(ax_onI, 'ON RF');
    colorbar(ax_onI);

    ax_onJ = subplot(5, 2, 8);
    imagesc(ax_onJ, onJ);
    colormap(ax_onJ, gray);
    axis(ax_onJ, 'off');
    title(ax_onJ, 'ON RF');
    colorbar(ax_onJ);

    % --- OFF RF ---
    ax_offI = subplot(5, 2, 9);
    imagesc(ax_offI, offI);
    colormap(ax_offI, gray);
    axis(ax_offI, 'off');
    title(ax_offI, 'OFF RF');
    colorbar(ax_offI);

    ax_offJ = subplot(5, 2, 10);
    imagesc(ax_offJ, offJ);
    colormap(ax_offJ, gray);
    axis(ax_offJ, 'off');
    title(ax_offJ, 'OFF RF');
    colorbar(ax_offJ);

    sgtitle(sprintf('Pair %d of %d  |  WAKE rank %d', figCount, nSummary, row), ...
        'FontSize', 11);

    % --- Save ---
    % savefig(fig, fullfile(saveDir, sprintf('pair_%02d_KS%d_KS%d.fig', figCount, ksIDi, ksIDj)));
    exportgraphics(fig, fullfile(saveDir, sprintf('pair_%02d_Cluster_%d_%d.png', figCount, ksIDi, ksIDj)), ...
        'Resolution', 150);
    close(fig);

    fprintf('Saved pair %d / %d (rank %d: KS%d -> KS%d)\n', figCount, nSummary, row, ksIDi, ksIDj);
end

fprintf('Done. %d summary figures saved to:\n  %s\n', figCount, saveDir);

%% --------------------------------------------------------------------------------
%% Compute overlap
%% --------------------------------------------------------------------------------

%% ---------- Compute RF similarity for all valid ranked pairs ----------

nPairs = height(ranked_wake_dedup);

% Preallocate result columns
corr_ON  = nan(nPairs, 1);
corr_OFF = nan(nPairs, 1);

jaccard_ON  = nan(nPairs, 1);
jaccard_OFF = nan(nPairs, 1);

centerDist_ON  = nan(nPairs, 1);
centerDist_OFF = nan(nPairs, 1);

for row = 1:nPairs
    ui = ranked_wake_dedup.ui(row);
    uj = ranked_wake_dedup.uj(row);

    % Skip invalid pairs
    if all(isnan(ccg_wake_norm(:,ui,uj))); continue; end

    % Get RF maps
    onI  = mean(RFmap{ui}.ON.OnSet,  3);
    offI = mean(RFmap{ui}.OFF.OnSet, 3);
    onJ  = mean(RFmap{uj}.ON.OnSet,  3);
    offJ = mean(RFmap{uj}.OFF.OnSet, 3);

    % --- 1. Pixel-wise correlation ---
    corr_ON(row)  = corr(onI(:),  onJ(:));
    corr_OFF(row) = corr(offI(:), offJ(:));

    % --- 2. Jaccard overlap of thresholded regions ---
    % Threshold: bins above 1 std of the map
    thresh_onI  = onI  > mean(onI(:))  + std(onI(:));
    thresh_onJ  = onJ  > mean(onJ(:))  + std(onJ(:));
    thresh_offI = offI > mean(offI(:)) + std(offI(:));
    thresh_offJ = offJ > mean(offJ(:)) + std(offJ(:));

    intersect_ON  = sum(thresh_onI(:)  & thresh_onJ(:));
    union_ON      = sum(thresh_onI(:)  | thresh_onJ(:));
    intersect_OFF = sum(thresh_offI(:) & thresh_offJ(:));
    union_OFF     = sum(thresh_offI(:) | thresh_offJ(:));

    if union_ON  > 0; jaccard_ON(row)  = intersect_ON  / union_ON;  end
    if union_OFF > 0; jaccard_OFF(row) = intersect_OFF / union_OFF; end

    % --- 3. Distance between RF centers (amplitude-weighted centroid) ---
    [nY, nX] = size(onI);
    [xGrid, yGrid] = meshgrid(1:nX, 1:nY);

    % ON centers
    wI = max(onI, 0);  wJ = max(onJ, 0);  % only positive weights
    if sum(wI(:)) > 0 && sum(wJ(:)) > 0
        cxI = sum(xGrid(:) .* wI(:)) / sum(wI(:));
        cyI = sum(yGrid(:) .* wI(:)) / sum(wI(:));
        cxJ = sum(xGrid(:) .* wJ(:)) / sum(wJ(:));
        cyJ = sum(yGrid(:) .* wJ(:)) / sum(wJ(:));
        centerDist_ON(row) = sqrt((cxI-cxJ)^2 + (cyI-cyJ)^2);  % pixels
    end

    % OFF centers
    wI = max(offI, 0);  wJ = max(offJ, 0);
    if sum(wI(:)) > 0 && sum(wJ(:)) > 0
        cxI = sum(xGrid(:) .* wI(:)) / sum(wI(:));
        cyI = sum(yGrid(:) .* wI(:)) / sum(wI(:));
        cxJ = sum(xGrid(:) .* wJ(:)) / sum(wJ(:));
        cyJ = sum(yGrid(:) .* wJ(:)) / sum(wJ(:));
        centerDist_OFF(row) = sqrt((cxI-cxJ)^2 + (cyI-cyJ)^2);  % pixels
    end
end

% Add to ranked table
ranked_wake_dedup.corr_ON       = corr_ON;
ranked_wake_dedup.corr_OFF      = corr_OFF;
ranked_wake_dedup.jaccard_ON    = jaccard_ON;
ranked_wake_dedup.jaccard_OFF   = jaccard_OFF;
ranked_wake_dedup.centerDist_ON  = centerDist_ON;
ranked_wake_dedup.centerDist_OFF = centerDist_OFF;

fprintf('RF similarity computed for %d pairs.\n', sum(~isnan(corr_ON)));

%% ---------- Plot RF similarity vs coupling strength ----------
figure('Units','centimeters','Position',[2 2 24 18]);

metrics    = {'corr_ON','corr_OFF','jaccard_ON','jaccard_OFF','centerDist_ON','centerDist_OFF'};
ylabels    = {'Pixel corr (ON)','Pixel corr (OFF)','Jaccard (ON)','Jaccard (OFF)','Center dist px (ON)','Center dist px (OFF)'};
color = [
    0.60 0.20 0.80;   % medium purple - corr_ON
    0.20 0.70 0.40;   % medium green  - corr_OFF
    0.40 0.10 0.60;   % dark purple   - jaccard_ON
    0.10 0.45 0.25;   % dark green    - jaccard_OFF
    0.80 0.60 0.95;   % light purple  - centerDist_ON
    0.65 0.90 0.65;   % light green   - centerDist_OFF
];

for m = 1:6
    ax = subplot(3, 2, m);
    validRows = ~isnan(ranked_wake_dedup.monoScore) & ~isnan(ranked_wake_dedup.(metrics{m}));

    xData = ranked_wake_dedup.(metrics{m})(validRows);
    yData = ranked_wake_dedup.monoScore(validRows);

    scatter(ax, xData, yData, 20, color(m,:), 'filled', 'MarkerFaceAlpha', 0.3);
    hold(ax, 'on');

    % Best fit line
    p    = polyfit(xData, yData, 1);
    xFit = linspace(min(xData), max(xData), 100);
    yFit = polyval(p, xFit);
    plot(ax, xFit, yFit, '-', 'Color', color(m,:), 'LineWidth', 2);

    xlabel(ax, ylabels{m});
    ylabel(ax, 'Coupling score (normalized)');
    r = corr(xData, yData, 'rows', 'complete');
    title(ax, sprintf('%s  |  r = %.2f', ylabels{m}, r));
    box(ax, 'off');
end
sgtitle('RF similarity vs WAKE coupling strength');

%% -------------- Plot receptive field overlap in summary figures ------------------



% %% Plot example pair with baseline-normalized CCG
% minBaseline  = 1.0;
% minRefSpikes = 50;
% 
% % Find first WAKE top pair valid across all states
% for row = 1:height(ranked_wake)
%     ui = ranked_wake.ui(row);
%     uj = ranked_wake.uj(row);
%     if ~all(isnan(ccg_nrem_rate(:,ui,uj))) && ~all(isnan(ccg_rem_rate(:,ui,uj)))
%         break
%     end
% end
% 
% [wake_norm, wake_ok, wake_reason] = normalize_ccg(ccg_wake_rate(:,ui,uj), flankIdx, minBaseline, minRefSpikes, nRefSpikes_wake(ui));
% [nrem_norm, nrem_ok, nrem_reason] = normalize_ccg(ccg_nrem_rate(:,ui,uj), flankIdx, minBaseline, minRefSpikes, nRefSpikes_nrem(ui));
% [rem_norm,  rem_ok,  rem_reason]  = normalize_ccg(ccg_rem_rate(:, ui,uj), flankIdx, minBaseline, minRefSpikes, nRefSpikes_rem(ui));
% 
% fprintf('WAKE: valid=%d  %s\n', wake_ok, wake_reason);
% fprintf('NREM: valid=%d  %s\n', nrem_ok, nrem_reason);
% fprintf('REM:  valid=%d  %s\n', rem_ok,  rem_reason);
% 
% figure; hold on;
% if wake_ok; plot(t_ccg*1000, wake_norm, 'k', 'LineWidth', 1.5); end
% if nrem_ok; plot(t_ccg*1000, nrem_norm, 'r', 'LineWidth', 1.5); end
% if rem_ok;  plot(t_ccg*1000, rem_norm,  'b', 'LineWidth', 1.5); end
% 
% yline(1, '--', 'Color', [0.5 0.5 0.5], 'Label', 'baseline');
% xline(0, '--', 'Color', [0.5 0.5 0.5]);
% xline([-5 5],  ':',  'Color', [0.7 0.7 0.7]);
% xlabel('Time lag (ms)');
% ylabel('Normalized rate (a.u.)');
% title(sprintf('CCG: KS%d → KS%d  (units %d, %d)', good_clusters(ui), good_clusters(uj), ui, uj));
% 
% legendEntries = {};
% if wake_ok; legendEntries{end+1} = 'WAKE'; end
% if nrem_ok; legendEntries{end+1} = 'NREM'; end
% if rem_ok;  legendEntries{end+1} = 'REM';  end
% legend(legendEntries, 'Location', 'best');
% box off;