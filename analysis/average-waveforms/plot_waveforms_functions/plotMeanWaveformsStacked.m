function plotMeanWaveformsStacked(meanWav, good_clusters, xpos, ypos, sr, clusterRange, ax)
% ax: axes handle to plot into

if nargin < 7 || isempty(ax)
    fig = figure; 
    ax = gca;
end

axes(ax);   % make sure we are plotting into the axes

[~, nsamp, ~] = size(meanWav);
num_channels_to_plot = 21;
spatialWindow = 300;
t = (1:nsamp) / sr * 1000;

for clusterID_index = clusterRange
    mean_wave_this_cluster = meanWav(:,:,clusterID_index);

    [~, peakCh] = max(max(abs(mean_wave_this_cluster),[],2));
    peakCh_ypos = ypos(peakCh);
    peakCh_xpos = xpos(peakCh);

    sameshank_channels = find(xpos == peakCh_xpos);
    peakParity = mod(peakCh, 2);
    sameParity_channels = sameshank_channels(mod(sameshank_channels,2)==peakParity);

    y_distance = abs(ypos(sameParity_channels) - peakCh_ypos);
    nearby_idx = y_distance <= spatialWindow;
    channelsToPlot = sameParity_channels(nearby_idx);

    if length(channelsToPlot) < num_channels_to_plot
        [~, sortIdx] = sort(y_distance);
        n_to_take = min(num_channels_to_plot, length(sameParity_channels));
        channelsToPlot = sameParity_channels(sortIdx(1:n_to_take));
    end
    if length(channelsToPlot) > num_channels_to_plot
        distances_to_plot = abs(ypos(channelsToPlot) - peakCh_ypos);
        [~, distIdx] = sort(distances_to_plot);
        channelsToPlot = channelsToPlot(distIdx(1:num_channels_to_plot));
    end
    [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
    channelsToPlot = channelsToPlot(sortIdx);

    offset = max(abs(mean_wave_this_cluster(:))) * 1.5;

    hold on;
    for i = 1:numel(channelsToPlot)
        ch = channelsToPlot(i);
        y_plot = mean_wave_this_cluster(ch,:) + (i-1)*offset;
        if ch == peakCh
            plot(t, y_plot, 'r', 'LineWidth', 2);
        else
            plot(t, y_plot, 'k', 'LineWidth', 1);
        end
    end

    % Add dashed horizontal line at peak channel's zero baseline
    peakIdxInPlot = find(channelsToPlot == peakCh);
    peakZeroY = (peakIdxInPlot - 1) * offset;
    yline(peakZeroY, ':', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);

    xlabel('Time (ms)');
    ylabel('Depth (superficial → deep)');
    title(sprintf('Unit %d | Peak Ch %d @ %.0f μm', ...
        good_clusters(clusterID_index), peakCh, peakCh_ypos));
    % yticks((0:numel(channelsToPlot)-1) * offset);
    yticklabels([]);
    % yticklabels(arrayfun(@(ch) sprintf('Ch %d (%.0f μm)', ch, ypos(ch)), ...
    %     channelsToPlot, 'UniformOutput', false));
end
