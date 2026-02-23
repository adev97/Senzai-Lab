function makePairSummaryPNGs(all_ccgs, all_ccgs_long, t, pairs, syn_pairs, good_clusters, ...
                              zJitter_wake, zJitter_nrem, zJitter_rem, ...
                              depthDiff, unitDepth, meanWav, cell_type, ...
                              xpos, ypos, sr, baselineBins_long, binSize_short, binSize_long, ...
                              pngDir)

% makePairSummaryPNGs - Create summary PNGs for jitter-confirmed synaptic pairs
%
% INPUTS:
%   all_ccgs          - cell array {ccg_wake, ccg_nrem, ccg_rem} (short CCGs)
%   all_ccgs_long     - cell array {ccg_wake_long, ccg_nrem_long, ccg_rem_long}
%   t                 - time vector from short CCG
%   pairs             - [nPairs x 2] matrix of ALL unit index pairs
%   syn_pairs         - indices into pairs to plot (e.g. syn_wake)
%   good_clusters     - original Kilosort cluster IDs
%   zJitter_wake/nrem/rem - jitter-corrected Z-scores for all pairs
%   depthDiff         - depth difference for each pair
%   unitDepth         - depth of each unit
%   meanWav           - mean waveforms [nChannels x nSamples x nClusters]
%   cell_type         - cell array of classification strings
%   xpos              - x positions of channels
%   ypos              - y positions of channels
%   sr                - sampling rate
%   baselineBins_long - baseline bin indices for long CCG normalization
%   binSize_short     - bin size of short CCG (s)
%   binSize_long      - bin size of long CCG (s)
%   pngDir            - output directory

if ~exist(pngDir, 'dir')
    mkdir(pngDir);
end

% Change orientation of cell_type if needed
if size(cell_type, 1) < size(cell_type, 2)
    cell_type = cell_type';
end

states = {'WAKE', 'NREM', 'REM'};
colors = {'k', 'r', 'b'};
nPairs = length(syn_pairs);

fprintf('Creating summary PNGs for %d synaptic pairs...\n', nPairs);

for j = 1:nPairs

    pairIdx  = syn_pairs(j);
    unitA    = pairs(pairIdx, 1);
    unitB    = pairs(pairIdx, 2);
    clusterA = good_clusters(unitA);
    clusterB = good_clusters(unitB);

    fprintf('  Pair %d/%d: Cluster %d & %d\n', j, nPairs, clusterA, clusterB);

    fig = figure('Position', [15, 100, 1200, 1275], 'Visible', 'off');

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

    % ===== Overall title =====
    sgtitle(sprintf('Pair #%d: Clusters %d & %d | Z_{wake}=%.1f | Z_{nrem}=%.1f | Z_{rem}=%.1f', ...
        j, clusterA, clusterB, ...
        zJitter_wake(pairIdx), zJitter_nrem(pairIdx), zJitter_rem(pairIdx)), ...
        'FontSize', 16, 'FontWeight', 'bold');

    % ===== Save PNG =====
    pngFilename = fullfile(pngDir, ...
        sprintf('Pair_%03d_Cluster_%d_%d.png', j, clusterA, clusterB));
    exportgraphics(fig, pngFilename, 'Resolution', 300);
    close(fig);

end

fprintf('\nDone! Created %d PNGs in:\n%s\n', nPairs, pngDir);
end