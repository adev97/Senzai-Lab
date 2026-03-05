function makePairSummaryWithRFsTHISONEWORKS(all_ccgs, t, pairs, good_clusters, ...
                                   coupling_wake, coupling_nrem, coupling_rem, ...
                                   unitDepth, meanWav, cell_type, ...
                                   xpos, ypos, sr, sortBy, topN, pdfDir, RFmap, maskVal, dispTime)

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

    %% ===== CCGs for all states =====
    subplot(7,2,[1 2]); hold on;
    plotMask = abs(t) <= dispTime;   % crop to ±dispTime
    for s = 1:3
        ccgPair = squeeze(all_ccgs{s}(:, unitA, unitB));
        plot(t(plotMask), ccgPair(plotMask), 'Color', colors{s}, 'LineWidth', 2);
    end
    xlabel('Time lag (s)');
    ylabel('Normalized CCG');
    title(sprintf('Cluster %d → Cluster %d', clusterA, clusterB));
    legend(states,'Location','best');
    grid on;
    xlim([-dispTime dispTime]);

    %% ===== Unit waveforms =====
    subplot(7,2,[3 5]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitA, gca);
    title(sprintf('Cluster %d', clusterA), 'FontSize',12,'FontWeight','bold');

    subplot(7,2,[4 6]);
    plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, unitB, gca);
    title(sprintf('Cluster %d', clusterB), 'FontSize',12,'FontWeight','bold');

    %% ===== Unit classification & depth =====
    subplot(7,2,7); axis off;
    text(0.5,0.7,sprintf('Cluster: %d', clusterA),'FontSize',10,'HorizontalAlignment','center');
    text(0.5,0.5,sprintf('Type: %s', cell_type{unitA}),'FontSize',11,'HorizontalAlignment','center','FontWeight','bold','Color','b');
    text(0.5,0.3,sprintf('Depth: %.1f μm', unitDepth(unitA)),'FontSize',10,'HorizontalAlignment','center');
    rectangle('Position',[0.05,0.05,0.9,0.9],'EdgeColor','k','LineWidth',1.5);

    subplot(7,2,8); axis off;
    text(0.5,0.7,sprintf('Cluster: %d', clusterB),'FontSize',10,'HorizontalAlignment','center');
    text(0.5,0.5,sprintf('Type: %s', cell_type{unitB}),'FontSize',11,'HorizontalAlignment','center','FontWeight','bold','Color','b');
    text(0.5,0.3,sprintf('Depth: %.1f μm', unitDepth(unitB)),'FontSize',10,'HorizontalAlignment','center');
    rectangle('Position',[0.05,0.05,0.9,0.9],'EdgeColor','k','LineWidth',1.5);

    %% ===== RF plots (unsmoothed) =====
    % Unit A
    subplot(7,2,9);
    rfA_on_raw  = mean(RFmap{unitA}.ON.OnSet,3) - RFmap{unitA}.baseline;
    imagesc(rfA_on_raw); axis image off; title(sprintf('Cluster %d - ON RF', clusterA)); colorbar

    subplot(7,2,11);
    rfA_off_raw = mean(RFmap{unitA}.OFF.OnSet,3) - RFmap{unitA}.baseline;
    imagesc(rfA_off_raw); axis image off; title(sprintf('Cluster %d - OFF RF', clusterA)); colorbar

    % Unit B
    subplot(7,2,10);
    rfB_on_raw  = mean(RFmap{unitB}.ON.OnSet,3) - RFmap{unitB}.baseline;
    imagesc(rfB_on_raw); axis image off; title(sprintf('Cluster %d - ON RF', clusterB)); colorbar

    subplot(7,2,12);
    rfB_off_raw = mean(RFmap{unitB}.OFF.OnSet,3) - RFmap{unitB}.baseline;
    imagesc(rfB_off_raw); axis image off; title(sprintf('Cluster %d - OFF RF', clusterB)); colorbar

    %% ===== Smoothed RFs for masks/overlaps only =====
    sigma = 1.2;
    hSize = 2 * ceil(2 * sigma) + 1;
    h = fspecial('gaussian', hSize, sigma);

    rfA_on_smooth  = imfilter(rfA_on_raw,  h, 'replicate');
    rfA_off_smooth = imfilter(rfA_off_raw, h, 'replicate');
    rfB_on_smooth  = imfilter(rfB_on_raw,  h, 'replicate');
    rfB_off_smooth = imfilter(rfB_off_raw, h, 'replicate');

    % Create logical masks
    noiseThresh_A_on  = 2 * std(rfA_on_smooth(:));
    noiseThresh_B_on  = 2 * std(rfB_on_smooth(:));
    noiseThresh_A_off = 2 * std(rfA_off_smooth(:));
    noiseThresh_B_off = 2 * std(rfB_off_smooth(:));

    mA_on  = (rfA_on_smooth  >= prctile(rfA_on_smooth(:),  maskVal)) & (rfA_on_smooth  > noiseThresh_A_on);
    mB_on  = (rfB_on_smooth  >= prctile(rfB_on_smooth(:),  maskVal)) & (rfB_on_smooth  > noiseThresh_B_on);
    mA_off = (rfA_off_smooth >= prctile(rfA_off_smooth(:), maskVal)) & (rfA_off_smooth > noiseThresh_A_off);
    mB_off = (rfB_off_smooth >= prctile(rfB_off_smooth(:), maskVal)) & (rfB_off_smooth > noiseThresh_B_off);

    %% ===== Overlap/Correlation plots =====
    subplot(7,2,13);
    onOverlay = zeros([size(mA_on),3]);
    onOverlay(:,:,1) = double(mA_on); 
    onOverlay(:,:,2) = double(mB_on); 
    imshow(onOverlay); axis image; title('ON Overlap (Yellow=Shared)');

    subplot(7,2,14);
    offOverlay = zeros([size(mA_off),3]);
    offOverlay(:,:,1) = double(mA_off); 
    offOverlay(:,:,2) = double(mB_off);
    imshow(offOverlay); axis image; title('OFF Overlap (Yellow=Shared)');

    %% ===== Title & Save =====
    sgtitle(sprintf('Pair #%d: Clusters %d & %d | %s Coupling: %.2f', ...
            j, clusterA, clusterB, sortLabel, couplingToSort(pairIdx)), 'FontSize',16,'FontWeight','bold');

    pngFilename = fullfile(pdfDir, sprintf('Pair_%03d_Cluster_%d_%d.png', j, clusterA, clusterB));
    exportgraphics(fig, pngFilename, 'Resolution',300);
    close(fig);
end

fprintf('Done! Created %d summary figures in %s\n', topN, pdfDir);
