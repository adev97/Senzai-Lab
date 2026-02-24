function makePairSummaryPNGsWORKS(all_ccgs, t, pairs, good_clusters, ...
                                   coupling_wake, coupling_nrem, coupling_rem, ...
                                   unitDepth, meanWav, cell_type, ...
                                   xpos, ypos, sr, sortBy, topN, pdfDir, RFmap)

% % Create summary PNGs for top N coupled pairs with RFs collapsed over time using mean
% makePairSummaryPNGsWORKS(all_ccgs, t, pairs, good_clusters, ...
%                         coupling_wake, coupling_nrem, coupling_rem, ...
%                         unitDepth, meanWav, cell_type, ...
%                         xpos, ypos, sr, sortBy, topN, pdfDir, RFmap);


% makePairSummaryPDFsSimple - Create summary PDFs/PNGs for top N coupled pairs
%
% INPUTS:
%   all_ccgs      - cell array {ccg_wake, ccg_nrem, ccg_rem} (already normalized)
%   t             - time vector from CCG
%   pairs         - [nPairs x 2] matrix of unit index pairs (pre-ranked)
%   good_clusters - original Kilosort cluster IDs
%   coupling_*    - normalized coupling values for each state
%   unitDepth     - depth of each unit
%   meanWav       - mean waveforms [nChannels x nSamples x nClusters]
%   cell_type     - cell array of classification strings
%   xpos, ypos    - x,y positions of channels
%   sr            - sampling rate
%   sortBy        - state to sort by: 'wake', 'nrem', 'rem'
%   topN          - number of top pairs to plot
%   pdfDir        - directory to save PDFs/PNGs

if nargin < 16
    error('Missing input arguments.');
end
if ~exist(pdfDir, 'dir')
    mkdir(pdfDir);
end

states = {'WAKE', 'NREM', 'REM'};
colors = {'k','r','b'};

% Select coupling vector for labeling
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

topN = min(topN, size(pairs,1));

for j = 1:topN
    pairIdx = j; % pairs are assumed pre-ranked
    unitA = pairs(pairIdx,1);
    unitB = pairs(pairIdx,2);
    clusterA = good_clusters(unitA);
    clusterB = good_clusters(unitB);

    fprintf('Creating summary for pair %d/%d: Cluster %d & %d...\n', ...
            j, topN, clusterA, clusterB);

    % Figure
    fig = figure('Position',[50,50,1200,1300]);

    % ===== CCGs for all states =====
    subplot(6,2,[1 2]); hold on;
    for s = 1:3
        ccgPair = squeeze(all_ccgs{s}(:, unitA, unitB));
        plot(t, ccgPair, 'Color', colors{s}, 'LineWidth', 2);
    end
    xlabel('Time lag (s)');
    ylabel('Normalized CCG');
    title(sprintf('Cluster %d → Cluster %d', clusterA, clusterB));
    legend(states,'Location','best');
    grid on;
    xlim([t(1) t(end)]);

    % ===== Unit A waveform =====
    subplot(6,2,[3 5]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitA, gca);
    title(sprintf('Cluster %d', clusterA), 'FontSize',12,'FontWeight','bold');

    % ===== Unit A classification & depth =====
    subplot(6,2,7); axis off;
    text(0.5,0.7,sprintf('Cluster: %d', clusterA),'FontSize',10,'HorizontalAlignment','center');
    text(0.5,0.5,sprintf('Type: %s', cell_type{unitA}),'FontSize',11,'HorizontalAlignment','center','FontWeight','bold','Color','b');
    text(0.5,0.3,sprintf('Depth: %.1f μm', unitDepth(unitA)),'FontSize',10,'HorizontalAlignment','center');
    rectangle('Position',[0.05,0.05,0.9,0.9],'EdgeColor','k','LineWidth',1.5);
   
    % ===== Unit A ON RF =====
    subplot(6,2,9);
    rfA_on = mean(RFmap{unitA}.ON.OnSet,3) - RFmap{unitA}.baseline;
    imagesc(rfA_on);
    axis image off
    title(sprintf('Cluster %d - ON RF', clusterA));
    colorbar
    
    % ===== Unit A OFF RF =====
    subplot(6,2,11);
    rfA_off = mean(RFmap{unitA}.OFF.OnSet,3) - RFmap{unitA}.baseline;
    imagesc(rfA_off);
    axis image off
    title(sprintf('Cluster %d -s OFF RF', clusterA));
    colorbar

    % ===== Unit B waveform =====
    subplot(6,2,[4 6]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitB, gca);
    title(sprintf('Cluster %d', clusterB), 'FontSize',12,'FontWeight','bold');

    % ===== Unit B classification & depth =====
    subplot(6,2,8); axis off;
    text(0.5,0.7,sprintf('Cluster: %d', clusterB),'FontSize',10,'HorizontalAlignment','center');
    text(0.5,0.5,sprintf('Type: %s', cell_type{unitB}),'FontSize',11,'HorizontalAlignment','center','FontWeight','bold','Color','b');
    text(0.5,0.3,sprintf('Depth: %.1f μm', unitDepth(unitB)),'FontSize',10,'HorizontalAlignment','center');
    rectangle('Position',[0.05,0.05,0.9,0.9],'EdgeColor','k','LineWidth',1.5);

    % ===== Unit B ON RF =====
    subplot(6,2,10);
    rfB_on = mean(RFmap{unitB}.ON.OnSet,3) - RFmap{unitB}.baseline;
    imagesc(rfB_on);
    axis image off
    title(sprintf('Cluster %d - ON RF', clusterB));
    colorbar

    % ===== Unit B OFF RF =====
    subplot(6,2,12);
    rfB_off = mean(RFmap{unitB}.OFF.OnSet,3) - RFmap{unitB}.baseline;
    imagesc(rfB_off);
    axis image off
    title(sprintf('Cluster %d - OFF RF', clusterB));
    colorbar

    % Overall title with coupling
    sgtitle(sprintf('Pair #%d: Clusters %d & %d | %s Coupling: %.2f', ...
            j, clusterA, clusterB, sortLabel, couplingToSort(pairIdx)), 'FontSize',16,'FontWeight','bold');

    % Save PNG
    pngFilename = fullfile(pdfDir, sprintf('Pair_%03d_Cluster_%d_%d.png', j, clusterA, clusterB));
    exportgraphics(fig, pngFilename, 'Resolution',300);

    % Optional: save PDF
    % pdfFilename = fullfile(pdfDir, sprintf('Pair_%03d_Cluster_%d_%d.pdf', j, clusterA, clusterB));
    % exportgraphics(fig, pdfFilename,'ContentType','vector','Resolution',300);

    close(fig);
end

fprintf('Done! Created %d summary figures in %s\n', topN, pdfDir);

end
