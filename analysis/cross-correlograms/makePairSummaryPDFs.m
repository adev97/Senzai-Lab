function makePairSummaryPDFs(all_ccgs, all_ccgs_long, t, pairs, good_clusters, ...
                              coupling_wake, coupling_nrem, coupling_rem, ...
                              depthDiff, unitDepth, meanWav, cell_type, ...
                              xpos, ypos, sr, baselineBins_long, binSize_short, binSize_long, ...
                              sortBy, topN, pdfDir)

% makePairSummaryPDFs - Create summary PDFs for top N coupled pairs
%
% INPUTS:
%   all_ccgs      - cell array {ccg_wake, ccg_nrem, ccg_rem} (short CCGs)
%   t             - time vector from short CCG
%   pairs         - [nPairs x 2] matrix of unit index pairs
%   good_clusters - original Kilosort cluster IDs
%   coupling_wake - normalized coupling strength during WAKE
%   coupling_nrem - normalized coupling strength during NREM
%   coupling_rem  - normalized coupling strength during REM
%   depthDiff     - depth difference for each pair
%   unitDepth     - depth of each unit
%   meanWav       - mean waveforms [nChannels x nSamples x nClusters]
%   cell_type     - cell array of classification strings (e.g. 'SC', 'RGC')
%   xpos          - x positions of channels
%   ypos          - y positions of channels
%   sr            - sampling rate
%   sortBy        - which state to sort by: 'wake', 'nrem', or 'rem'
%   topN          - how many pairs to create PDFs for
%   pdfDir        - directory to save PDFs

% Input validation
if nargin < 21
    error('Missing input arguments. Check function signature.');
end

if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end

% Select coupling vector to sort by
switch lower(sortBy)
    case 'wake'
        couplingToSort = coupling_wake;
        sortLabel = 'WAKE';
    case 'nrem'
        couplingToSort = coupling_nrem;
        sortLabel = 'NREM';
    case 'rem'
        couplingToSort = coupling_rem;
        sortLabel = 'REM';
    otherwise
        error('sortBy must be ''wake'', ''nrem'', or ''rem''');
end

% Sort pairs by chosen coupling
[~, sortIdx] = sort(couplingToSort, 'descend');
pairsToPlot  = sortIdx(1:topN);

% Normalization parameters for visualization
% nBins     = size(all_ccgs{1}, 1);
% nEdgeBins = floor(0.2 * nBins);
states    = {'WAKE', 'NREM', 'REM'};
colors    = {'k', 'r', 'b'};

% Change orientation of cell_type if needed
if size(cell_type, 1) < size(cell_type, 2)
    cell_type = cell_type';
end

%% Loop through pairs and create PDF
for j = 1:length(pairsToPlot)
    pairIdx = pairsToPlot(j);

    unitA    = pairs(pairIdx, 1);
    unitB    = pairs(pairIdx, 2);
    clusterA = good_clusters(unitA);
    clusterB = good_clusters(unitB);

    fprintf('Creating PDF for pair %d/%d: Cluster %d & %d...\n', ...
        j, length(pairsToPlot), clusterA, clusterB);

    % Create figure
    fig = figure('Position', [15, 100, 1200, 1275]);

    % ===== ROW 1: CCG across all states =====
    subplot(6, 2, [1 2])
    hold on;

    for s = 1:length(states)
        ccgState = all_ccgs{s};
        ccgPair  = squeeze(ccgState(:, unitA, unitB));
        ccgLong         = squeeze(all_ccgs_long{s}(:, unitA, unitB));
        baseline        = median(ccgLong(baselineBins_long));
        baseline_scaled = baseline * (binSize_short / binSize_long);
        if baseline_scaled == 0
            ccgNorm = ccgPair;
        else
            ccgNorm = ccgPair / baseline_scaled;
        end

        % Normalize using edge bins of short CCG
        % baselineBins = [1:nEdgeBins, nBins-nEdgeBins+1:nBins];
        % baseline     = median(ccgPair(baselineBins));
        % 
        % if baseline == 0
        %     ccgNorm = ccgPair;
        % else
        %     ccgNorm = ccgPair / baseline;
        % end

        plot(t, ccgNorm, 'Color', colors{s}, 'LineWidth', 2);
    end

    xlabel('Time lag (s)');
    ylabel('Normalized CCG');
    title(sprintf('Cluster %d → Cluster %d (Units %d → %d)', ...
        clusterA, clusterB, unitA, unitB));
    legend(states, 'Location', 'best');
    grid on;
    xlim([t(1) t(end)]);

    % ===== ROW 2-3: Unit A Waveform =====
    subplot(6, 2, [3 5]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitA, gca);
    title(sprintf('Cluster %d', clusterA), 'FontSize', 12, 'FontWeight', 'bold');

    % ===== ROW 4: Unit A Classification =====
    subplot(6, 2, 7)
    axis off;
    text(0.5, 0.7, sprintf('Cluster: %d', clusterA), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    text(0.5, 0.5, sprintf('Type: %s', cell_type{unitA}), 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
        'Color', 'b', 'VerticalAlignment', 'middle');
    text(0.5, 0.3, sprintf('Depth: %.1f μm', unitDepth(unitA)), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.05, 0.05, 0.9, 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);

    % ===== ROW 5-6: Unit A Receptive Field =====
    subplot(6, 2, [9 11])
    axis off;
    text(0.5, 0.5, sprintf('Receptive Field\nUnit A\n\n[Add RF data]'), ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.1, 0.1, 0.8, 0.8], 'EdgeColor', 'k', 'LineStyle', '--');

    % ===== ROW 2-3: Unit B Waveform =====
    subplot(6, 2, [4 6]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitB, gca);
    title(sprintf('Cluster %d', clusterB), 'FontSize', 12, 'FontWeight', 'bold');

    % ===== ROW 4: Unit B Classification =====
    subplot(6, 2, 8)
    axis off;
    text(0.5, 0.7, sprintf('Cluster: %d', clusterB), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    text(0.5, 0.5, sprintf('Type: %s', cell_type{unitB}), 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
        'Color', 'b', 'VerticalAlignment', 'middle');
    text(0.5, 0.3, sprintf('Depth: %.1f μm', unitDepth(unitB)), 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.05, 0.05, 0.9, 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);

    % ===== ROW 5-6: Unit B Receptive Field =====
    subplot(6, 2, [10 12])
    axis off;
    text(0.5, 0.5, sprintf('Receptive Field\nUnit B\n\n[Add RF data]'), ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    rectangle('Position', [0.1, 0.1, 0.8, 0.8], 'EdgeColor', 'k', 'LineStyle', '--');

    % Overall title
    sgtitle(sprintf('Pair #%d: Clusters %d & %d | %s Coupling: %.2f', ...
        j, clusterA, clusterB, sortLabel, couplingToSort(pairIdx)), ...
        'FontSize', 16, 'FontWeight', 'bold');

    % % Save PDF
    % pdfFilename = fullfile(pdfDir, ...
    %     sprintf('Pair_%03d_Cluster_%d_%d.pdf', j, clusterA, clusterB));
    % exportgraphics(fig, pdfFilename, 'ContentType', 'vector', 'Resolution', 300);

    % Save PNG
    pngFilename = fullfile(pdfDir, ...
    sprintf('Pair_%03d_Cluster_%d_%d.png', j, clusterA, clusterB));

    exportgraphics(fig, pngFilename, 'Resolution', 300);

    close(fig);
end

fprintf('\nDone! Created %d PDFs in:\n%s\n', topN, pdfDir);
end