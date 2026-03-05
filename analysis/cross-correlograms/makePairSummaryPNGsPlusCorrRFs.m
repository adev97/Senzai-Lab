function makePairSummaryPNGsPlusCorrRFs(all_ccgs, t, pairs, good_clusters, ...
                                   coupling_wake, coupling_nrem, coupling_rem, ...
                                   unitDepth, meanWav, cell_type, ...
                                   xpos, ypos, sr, sortBy, topN, pdfDir, RFmap, maskVal, dispTime)

% % Create summary PNGs for top N coupled pairs with RFs collapsed over time using mean
% makePairSummaryPNGsWORKS(all_ccgs, t, pairs, good_clusters, ...
%                         coupling_wake, coupling_nrem, coupling_rem, ...
%                         unitDepth, meanWav, cell_type, ...
%                         xpos, ypos, sr, sortBy, topN, pdfDir, RFmap);

% Uses MEAN to collapse the RFs


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
    fig = figure('Position',[50,50,1200,1500], 'Visible','off');

    % ===== CCGs for all states =====
    subplot(7,2,[1 2]); hold on;
    
    % Define plot mask for ±0.5 s
    plotMask = abs(t) <= dispTime;   % crop to ±0.5 s
    
    for s = 1:3
        ccgPair = squeeze(all_ccgs{s}(:, unitA, unitB));
        plot(t(plotMask), ccgPair(plotMask), 'Color', colors{s}, 'LineWidth', 2);
    end
    
    xlabel('Time lag (s)');
    ylabel('Normalized CCG');
    title(sprintf('Cluster %d → Cluster %d', clusterA, clusterB));
    legend(states,'Location','best');
    grid on;
    xlim([-dispTime dispTime]);   % crop axes to ±0.5 s

    % ===== Unit A waveform =====
    subplot(7,2,[3 5]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitA, gca);
    title(sprintf('Cluster %d', clusterA), 'FontSize',12,'FontWeight','bold');

    % ===== Unit A classification & depth =====
    subplot(7,2,7); axis off;
    text(0.5,0.7,sprintf('Cluster: %d', clusterA),'FontSize',10,'HorizontalAlignment','center');
    text(0.5,0.5,sprintf('Type: %s', cell_type{unitA}),'FontSize',11,'HorizontalAlignment','center','FontWeight','bold','Color','b');
    text(0.5,0.3,sprintf('Depth: %.1f μm', unitDepth(unitA)),'FontSize',10,'HorizontalAlignment','center');
    rectangle('Position',[0.05,0.05,0.9,0.9],'EdgeColor','k','LineWidth',1.5);
   
    % ===== Unit A ON RF =====
    subplot(7,2,9);
    rfA_on = mean(RFmap{unitA}.ON.OnSet,3) - RFmap{unitA}.baseline;
    imagesc(rfA_on);
    axis image off
    title(sprintf('Cluster %d - ON RF', clusterA));
    colorbar
    
    % ===== Unit A OFF RF =====
    subplot(7,2,11);
    rfA_off = mean(RFmap{unitA}.OFF.OnSet,3) - RFmap{unitA}.baseline;
    imagesc(rfA_off);
    axis image off
    title(sprintf('Cluster %d - OFF RF', clusterA));
    colorbar

    % ===== Unit B waveform =====
    subplot(7,2,[4 6]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitB, gca);
    title(sprintf('Cluster %d', clusterB), 'FontSize',12,'FontWeight','bold');

    % ===== Unit B classification & depth =====
    subplot(7,2,8); axis off;
    text(0.5,0.7,sprintf('Cluster: %d', clusterB),'FontSize',10,'HorizontalAlignment','center');
    text(0.5,0.5,sprintf('Type: %s', cell_type{unitB}),'FontSize',11,'HorizontalAlignment','center','FontWeight','bold','Color','b');
    text(0.5,0.3,sprintf('Depth: %.1f μm', unitDepth(unitB)),'FontSize',10,'HorizontalAlignment','center');
    rectangle('Position',[0.05,0.05,0.9,0.9],'EdgeColor','k','LineWidth',1.5);

    % ===== Unit B ON RF =====
    subplot(7,2,10);
    rfB_on = mean(RFmap{unitB}.ON.OnSet,3) - RFmap{unitB}.baseline;
    imagesc(rfB_on);
    axis image off
    title(sprintf('Cluster %d - ON RF', clusterB));
    colorbar

    % ===== Unit B OFF RF =====
    subplot(7,2,12);
    rfB_off = mean(RFmap{unitB}.OFF.OnSet,3) - RFmap{unitB}.baseline;
    imagesc(rfB_off);
    axis image off
    title(sprintf('Cluster %d - OFF RF', clusterB));
    colorbar

    % --- Inside the loop, before calculating masks ---
    sigma = 1.2;
    % Create a Gaussian filter kernel (Bypasses Control System Toolbox requirement)
    hSize = 2 * ceil(2 * sigma) + 1; % Standard window size for sigma
    h = fspecial('gaussian', hSize, sigma); 

    % 1. Process Unit A
    rawA_on  = mean(RFmap{unitA}.ON.OnSet, 3)  - RFmap{unitA}.baseline;
    rawA_off = mean(RFmap{unitA}.OFF.OnSet, 3) - RFmap{unitA}.baseline;
    rfA_on   = imfilter(rawA_on, h, 'replicate');
    rfA_off  = imfilter(rawA_off, h, 'replicate');
    
    % 2. Process Unit B
    rawB_on  = mean(RFmap{unitB}.ON.OnSet, 3)  - RFmap{unitB}.baseline;
    rawB_off = mean(RFmap{unitB}.OFF.OnSet, 3) - RFmap{unitB}.baseline;
    rfB_on   = imfilter(rawB_on, h, 'replicate');
    rfB_off  = imfilter(rawB_off, h, 'replicate');
    
    % 3. Generate Logical Masks (Top 20%)
    % mA_on = rfA_on >= prctile(rfA_on(:), maskVal);
    % mB_on = rfB_on >= prctile(rfB_on(:), maskVal);
    % mA_off = rfA_off >= prctile(rfA_off(:), maskVal);
    % mB_off = rfB_off >= prctile(rfB_off(:), maskVal);
    noiseThresh_A_on  = 2 * std(rfA_on(:));
    noiseThresh_B_on  = 2 * std(rfB_on(:));
    noiseThresh_A_off = 2 * std(rfA_off(:));
    noiseThresh_B_off = 2 * std(rfB_off(:));
    
    mA_on  = (rfA_on  >= prctile(rfA_on(:),  maskVal)) & (rfA_on  > noiseThresh_A_on);
    mB_on  = (rfB_on  >= prctile(rfB_on(:),  maskVal)) & (rfB_on  > noiseThresh_B_on);
    mA_off = (rfA_off >= prctile(rfA_off(:), maskVal)) & (rfA_off > noiseThresh_A_off);
    mB_off = (rfB_off >= prctile(rfB_off(:), maskVal)) & (rfB_off > noiseThresh_B_off);

    % --- Row 7: The Overlap/Correlation Row ---
    % Subplot 13: ON Overlap (Red=A, Green=B, Yellow=Overlap)
    subplot(7,2,13);
    onOverlay = zeros([size(mA_on), 3]);
    onOverlay(:,:,1) = double(mA_on); % Convert to double for imshow compatibility
    onOverlay(:,:,2) = double(mB_on); 
    imshow(onOverlay);
    axis image; title('ON Overlap (Yellow=Shared)');
    
    % Subplot 14: OFF Overlap
    subplot(7,2,14);
    offOverlay = zeros([size(mA_off), 3]);
    offOverlay(:,:,1) = double(mA_off); 
    offOverlay(:,:,2) = double(mB_off);
    imshow(offOverlay);
    axis image; title('OFF Overlap (Yellow=Shared)');

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
