%% Script to build summary pdfs that contain each pair's ccgs (50ms lag), each unit's average waveform, and each unit's classification

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master")) 

sr = 30000;
nchan_probe = 384;
dtype = 'int16';

% Constants
nClusters = length(good_clusters);

% kilosort directory
ksDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Elissa_Belluccini\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';

spikeTimes    = readNPY(fullfile(ksDir,'spike_times.npy'));
spikeTimes    = spikeTimes + 1;
spikeClusters = readNPY(fullfile(ksDir,'spike_clusters.npy')); % loads all units (not only the good ones)
templates     = readNPY(fullfile(ksDir,'templates.npy')); % can use template to compare? see if average waveforms from raw data are super different from  template
chanPos       = readNPY(fullfile(ksDir,'channel_positions.npy')); % [chan × 2]

% Directory to save PDFs
pdfDir = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\PAIR_SUMMARIES";
if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end

%% Load Mean Waveform Data
meanWaveformDir = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters";
load(fullfile(meanWaveformDir, 'meanWav_units.mat'));

%% Load RGC-SC Classification
classificationDir = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\SC-RGC-classification-with-kmeans-all-channels-BEST";
load(fullfile(classificationDir, 'RGC-SC-classification.mat'));

%% Load CCGs
ccgDir = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\ccgs-all-units-40ms-Lag";
load(fullfile(ccgDir, 'ccgs_all.mat'));

% get unit depth
unitDepth = zeros(nClusters,1);

for i = 1:nClusters
    clu = good_clusters(i);
    tempIdx = mode(spikeClusters(spikeClusters==clu));
    template = squeeze(templates(tempIdx+1,:,:)); % ks indexing
    
    [~,peakChan] = max(max(abs(template),[],2));
    unitDepth(i) = ypos(peakChan);
end

% Do coupling calculation, build pairs
pairs = nchoosek(1:nClusters,2);
nPairs = size(pairs,1);
depthDiff = zeros(nPairs,1);

% extract coupling strength per pair of units based on depth - same as make_ccgs
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

states = {'WAKE', 'NREM', 'REM'};
all_ccgs = {ccg_wake, ccg_nrem, ccg_rem};
nStates = length(states);

nBins = size(ccg_wake, 1);

% change orientation for cell type classification arrays
cell_type = cell_type';
cluster_id = cluster_id';

% Select which pairs to create summaries for
topN = 10; % create PDFs for top 10 pairs
[~, idx_wake] = sort(coupling_wake, 'descend'); % from WAKE, top 10 most coupled
pairsToPlot = idx_wake(1:topN);

%% Loop through pairs and create PDF

for i = 1:length(pairsToPlot)
    pairIdx = pairsToPlot(i);

    unitA = pairs(pairIdx, 1);
    unitB = pairs(pairIdx, 2);
    clusterA = good_clusters(unitA);
    clusterB = good_clusters(unitB);

    meanWaveformA = meanWav(unitA, :);
    meanWaveformB = meanWav(unitB, :);

    fprintf('Creating PDF for pair %d/%d: Cluster %d & %d...\n', i, length(pairsToPlot), clusterA, clusterB);

    % Create a new PDF for the current pair
    % pdfFileName = fullfile(pdfDir, sprintf('Pair_%d.pdf', i));
    % pdf = pdfwrite(pdfFileName);

    % Create figure
    fig = figure('Position', [15, 100, 1200, 1275]);

    % ------ Panel 1: CCGs across all states (wake, nrem, rem)
    subplot(6, 2, [1 2])
    hold on;

    % Check if clusters exist in your good units
    if isempty(unitA) || isempty(unitB)
        error('One or both cluster IDs not found in good_clusters');
    end
    
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

    % ----- Panel 2: Unit A Waveform
    subplot(6, 2, [3 5]);
    % hold on;
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitA, gca);
    title(sprintf('Cluster %d', clusterA), 'FontSize', 12, 'FontWeight', 'bold');


    % --- Unit A classification metrics
    subplot(6, 2, 7)
    axis off;
    text(0.5, 0.7, sprintf('Cluster: %d', clusterA), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    text(0.5, 0.5, sprintf('Type: %s', cell_type{unitA}), 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'b', 'VerticalAlignment', 'middle');
    text(0.5, 0.3, sprintf('Depth: %.1f μm', unitDepth(unitA)), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.05, 0.05, 0.9, 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    % Unit A Receptive Field
    subplot(6, 2, [9 11])
    axis off;
    text(0.5, 0.5, sprintf('Receptive Field\nUnit A\n\n[Add RF data]'), ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.1, 0.1, 0.8, 0.8], 'EdgeColor', 'k', 'LineStyle', '--');
    

    %  % ===== MIDDLE COLUMN: Coupling Metrics =====
    % 
    % subplot(4, 3, [5 8 11])
    % axis off;
    % 
    % textStr = {
    %     sprintf('Coupling Metrics'),
    %     '',
    %     sprintf('WAKE:  %.1f', coupling_wake(pairIdx)),
    %     sprintf('NREM:  %.1f', coupling_nrem(pairIdx)),
    %     sprintf('REM:   %.1f', coupling_rem(pairIdx)),
    %     '',
    %     sprintf('NREM - WAKE: %.1f', coupling_nrem(pairIdx) - coupling_wake(pairIdx)),
    %     sprintf('REM - WAKE:  %.1f', coupling_rem(pairIdx) - coupling_wake(pairIdx)),
    %     '',
    %     sprintf('Depth diff: %.1f μm', depthDiff(pairIdx)),
    %     '',
    %     sprintf('Spike counts:'),
    %     sprintf('  Unit A: %d', sum(spikeClusters == clusterA)),
    %     sprintf('  Unit B: %d', sum(spikeClusters == clusterB))
    % };
    % 
    % text(0.1, 0.95, textStr, 'FontSize', 11, 'VerticalAlignment', 'top', ...
    %     'FontName', 'FixedWidth', 'FontWeight', 'bold');
    % rectangle('Position', [0.05, 0.05, 0.9, 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);
    

    % ----- Panel 3: Unit B Waveform
    subplot(6, 2, [4 6]);
    % hold on;
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitB, gca);
    title(sprintf('Cluster %d', clusterB), 'FontSize', 12, 'FontWeight', 'bold');

    % Unit B classification
    subplot(6, 2, 8)
    axis off;
    text(0.5, 0.7, sprintf('Cluster: %d', clusterB), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    text(0.5, 0.5, sprintf('Type: %s', cell_type{unitB}), 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'b', 'VerticalAlignment', 'middle');
    text(0.5, 0.3, sprintf('Depth: %.1f μm', unitDepth(unitB)), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.05, 0.05, 0.9, 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);   
    
    % Unit B Receptive Field
    subplot(6, 2, [10 12])
    axis off;
    text(0.5, 0.5, sprintf('Receptive Field\nUnit B\n\n[Add RF data]'), ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.1, 0.1, 0.8, 0.8], 'EdgeColor', 'k', 'LineStyle', '--');
    
    % Add overall title
    sgtitle(sprintf('Pair Summary #%d: Clusters %d & %d | WAKE Coupling: %.1f', ...
        p, clusterA, clusterB, coupling_nrem(pairIdx)), ...
        'FontSize', 16, 'FontWeight', 'bold');
    
end




