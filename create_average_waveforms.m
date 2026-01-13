%% New Script that modifies DepthSort_meanWaveForms (from Yuta) for current Neuropixel data and only makes average waveform, without sorting

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("D:\buzcode-master")) 

sr = 30000;
nchan_probe = 384;
dtype = 'int16';

ksDir = 'D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';

spikeTimes    = readNPY(fullfile(ksDir,'spike_times.npy'));
spikeTimes    = spikeTimes + 1;
spikeClusters = readNPY(fullfile(ksDir,'spike_clusters.npy')); % loads all units (not only the good ones)
templates     = readNPY(fullfile(ksDir,'templates.npy')); % can use template to compare? see if average waveforms from raw data are super different from  template
chanPos       = readNPY(fullfile(ksDir,'channel_positions.npy')); % [chan × 2]

% load cluster description (good/mua/noise)
cgFile = fullfile(ksDir,'cluster_group.tsv');

cluster_groups = readtable(cgFile, ...
    'FileType','text', ...
    'Delimiter','\t');

keepGroups = {'good'}; % only keep good labeled units

toKeep = ismember(cluster_groups.group, keepGroups);
keepClusters = cluster_groups.cluster_id(toKeep); % contains only the good cluster ids (162 total)

% filter spikes by only kept clusters
keepSpike = ismember(spikeClusters, keepClusters);

spikeTimes    = spikeTimes(keepSpike);
spikeClusters = spikeClusters(keepSpike);

% fill x position and y position per channel
xpos = chanPos(:,1);
ypos = chanPos(:,2);

good_clusters = unique(spikeClusters);
good_clusters(good_clusters==0) = [];   % remove noise if present
nClusters = numel(good_clusters);

% Begin waveform extraction

rawFile = "D:\Kilosort\Mouse08_SC_20251007_810to2250\Merged.dat";


% chunking it because Merged.dat is huge
nCh = 385; % from Elissa
bytesPerSample = 2;

fileInfo = dir(rawFile);
nSamples = fileInfo.bytes / (bytesPerSample * nCh);

sbefore = 30;   % samples before trough (-1 ms)
safter  = 60;   % samples after trough (+2 ms)
nsamp   = sbefore + safter;

maxSpikes = 2000;   % to sample per cluster

% pre-allocate
meanWav = zeros(nchan_probe, nsamp, nClusters);

fid = fopen(rawFile,'r');

for c = 1:nClusters
    cluID = good_clusters(c);
    idx = find(spikeClusters == cluID);

    if numel(idx) > maxSpikes
        idx = idx(round(linspace(1,numel(idx),maxSpikes)));
    end

    wav = zeros(384, nsamp, numel(idx));

    for s = 1:numel(idx)

        t = spikeTimes(idx(s));

        if t <= sbefore || t + safter > nSamples
            continue
        end

        % Determine which chunk the spike is in
        fseek(fid, (t - sbefore - 1) * nCh * bytesPerSample, 'bof');
        raw = fread(fid, nCh * nsamp, 'int16=>double');
        w = reshape(raw, [nCh, nsamp]);
        w = w(1:nchan_probe, :);
        % w = fread(fid, [nCh, nsamp], '*int16');
        w = w - median(w(:,1:sbefore),2); % baseline
        wav(:,:,s) = w;
    end

    meanWav(:,:,c) = mean(wav,3);
end

fclose(fid);

disp("Calculate Avg Waveforms Complete");

% now you have the mean wave form for each good unit!

%% Visualize mean waveforms - filter by channel parity to avoid checkerboard (even vs odd units)

outWaveDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\plus2ms-minus1ms\average-waveforms-per-unit-from-raw';
outDepthDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\plus2ms-minus1ms\heatmap-waveforms-per-unit-from-raw';

if ~exist(outWaveDir, 'dir'); mkdir(outWaveDir); end
if ~exist(outDepthDir, 'dir'); mkdir(outDepthDir); end


for clusterID_index = 1:nClusters

    mean_wave_this_cluster = meanWav(:,:,clusterID_index);
    
    % Find channel with largest spike amplitude
    [~, peakCh] = max(max(abs(mean_wave_this_cluster),[],2));
    
    peakCh_ypos = ypos(peakCh);
    peakCh_xpos = xpos(peakCh);
    
    % Get channels on same shank
    sameShank = (xpos == peakCh_xpos);
    
    % Filter by same parity (odd/even) as peak channel
    peakParity = mod(peakCh, 2);  % 0 for even, 1 for odd
    sameParity = (mod(find(sameShank), 2) == peakParity);
    
    % Get indices in the original channel array
    sameshank_channels = find(sameShank);
    sameParity_channels = sameshank_channels(sameParity);
    
    % Calculate distance from peak channel (in y-direction)
    y_distance = abs(ypos(sameParity_channels) - peakCh_ypos);
    
    % Start with distance window
    spatialWindow = 150; % micrometers
    nearby_idx = find(y_distance <= spatialWindow);
    channelsToPlot = sameParity_channels(nearby_idx);
    
    % Sort by depth
    [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
    channelsToPlot = channelsToPlot(sortIdx);
    
    % If too few channels, expand to nearest neighbors with same parity
    if length(channelsToPlot) < 11
        warning('Only %d channels within %d μm with same parity, expanding to nearest neighbors', ...
                length(channelsToPlot), spatialWindow);
        
        [~, sortIdx] = sort(y_distance);
        sorted_channels = sameParity_channels(sortIdx);
        
        % Take up to 11 nearest channels
        n_to_take = min(11, length(sorted_channels));
        channelsToPlot = sorted_channels(1:n_to_take);
        
        [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
        channelsToPlot = channelsToPlot(sortIdx);
    end
    
    % If too many, keep closest to peak
    if length(channelsToPlot) > 11
        distances_to_plot = abs(ypos(channelsToPlot) - peakCh_ypos);
        [~, distIdx] = sort(distances_to_plot);
        channelsToPlot = channelsToPlot(distIdx(1:11));
        
        [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
        channelsToPlot = channelsToPlot(sortIdx);
    end
    
    % Plot waveforms stacked by depth
    fig1 = figure('Visible','off');
    hold on;
    
    t = (1:nsamp) / sr * 1000;   % time in ms
    offset = max(abs(mean_wave_this_cluster(:))) * 1.5;  % vertical spacing
    
    for i = 1:numel(channelsToPlot)
        ch = channelsToPlot(i);
        y_plot = mean_wave_this_cluster(ch,:) + (i-1)*offset;
        
        if ch == peakCh
            plot(t, y_plot, 'r', 'LineWidth', 2.); % peak channel in red
        else
            plot(t, y_plot, 'k', 'LineWidth', 1);
        end
    end
    
    % Add formatting
    xlabel('Time (ms)');
    ylabel('Depth (superficial → deep)');
    parityStr = {'even', 'odd'};
    title(sprintf('Unit %d | Peak Ch %d @ %.0f μm (%s)', ...
            good_clusters(clusterID_index), peakCh, peakCh_ypos, ...
            parityStr{peakParity+1}));
     
    peakIdxInPlot = find(channelsToPlot == peakCh);
    peakZeroY = (peakIdxInPlot - 1) * offset;
    yline(peakZeroY, ':', 'LineWidth', 0.75, 'Color', [0.4 0.4 0.4]);
    
    
    % Label with actual channel numbers and depths
    yticks((0:numel(channelsToPlot)-1) * offset);
    yticklabels(arrayfun(@(ch) sprintf('Ch %d (%.0f μm)', ch, ypos(ch)), ...
                channelsToPlot, 'UniformOutput', false));
    
    % % Add text annotation
    % text(0.02, 0.98, sprintf('Shank x=%.0f μm, %s channels\nSpatial window: ±%.0f μm', ...
    %      peakCh_xpos, parityStr{peakParity+1}, spatialWindow), ...
    %      'Units', 'normalized', 'VerticalAlignment', 'top', ...
    %      'FontSize', 10, 'BackgroundColor', 'white');
    saveas(fig1, fullfile(outWaveDir, ...
        sprintf('unit_%03d_avg.png', good_clusters(clusterID_index))));
    close(fig1);

    hold off;
    
    
    %% Check depth sorting by visualizing spike propagation
    
    fig2 = figure('Visible','off');
    subplot(1,2,1);
    imagesc(t, 1:length(channelsToPlot), mean_wave_this_cluster(channelsToPlot,:));
    colormap('jet');
    colorbar;
    xlabel('Time (ms)');
    ylabel('Channel (superficial → deep)');
    title('Waveform heatmap');
    set(gca, 'YTick', 1:length(channelsToPlot));
    set(gca, 'YTickLabel', arrayfun(@(ch) sprintf('Ch%d', ch), channelsToPlot, 'UniformOutput', false));
    
    % Find peak time on each channel
    subplot(1,2,2);
    peak_times = zeros(length(channelsToPlot), 1);
    for i = 1:length(channelsToPlot)
        ch = channelsToPlot(i);
        [~, peak_idx] = min(mean_wave_this_cluster(ch,:)); % trough time
        peak_times(i) = t(peak_idx);
    end
    
    plot(peak_times, ypos(channelsToPlot), 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    xlabel('Trough time (ms)');
    ylabel('Depth (μm)');
    title('Spike propagation');
    grid on;
    
    saveas(fig2, fullfile(outDepthDir, ...
        sprintf('unit_%03d_heatmap.png', good_clusters(clusterID_index))));
    close(fig2);
    
    % Should see a smooth progression if sorting is correct
end







































% %% Visualize mean waveforms with spatial information preserved
% 
% clusterID_index = 103;  % which cluster INDEX (not cluster ID!)
% mean_wave_this_cluster = meanWav(:,:,clusterID_index);
% 
% % Find channel with largest spike amplitude
% [~, peakCh] = max(max(abs(mean_wave_this_cluster),[],2));
% 
% % Get the spatial position of the peak channel
% peakCh_ypos = ypos(peakCh);
% peakCh_xpos = xpos(peakCh);
% 
% % Get channels on the same shank AND same column (left/right)
% sameShank = (xpos == peakCh_xpos);  % This already filters by exact x position!
% 
% % Calculate distance from peak channel (in y-direction ONLY, since x is already matched)
% y_distance = abs(ypos - peakCh_ypos);
% 
% % For same-column channels only
% samecol_channels = find(sameShank);
% 
% % Start with distance window
% spatialWindow = 150; % micrometers - adjust based on your needs
% nearby = sameShank & (y_distance <= spatialWindow);
% channelsToPlot = find(nearby);
% 
% % Sort by depth
% [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
% channelsToPlot = channelsToPlot(sortIdx);
% 
% % If too few channels, expand to nearest neighbors ON SAME COLUMN
% if length(channelsToPlot) < 11
%     warning('Only %d channels within %d μm on same column, expanding to nearest neighbors', ...
%             length(channelsToPlot), spatialWindow);
% 
%     [~, sortIdx] = sort(y_distance(samecol_channels));
%     sorted_channels = samecol_channels(sortIdx);
% 
%     % Take up to 11 nearest channels
%     n_to_take = min(11, length(sorted_channels));
%     channelsToPlot = sorted_channels(1:n_to_take);
% 
%     [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
%     channelsToPlot = channelsToPlot(sortIdx);
% end
% 
% % If too many, keep closest to peak
% if length(channelsToPlot) > 11
%     [~, distIdx] = sort(y_distance(channelsToPlot));
%     channelsToPlot = channelsToPlot(distIdx(1:11));
% 
%     [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
%     channelsToPlot = channelsToPlot(sortIdx);
% end
% 
% % Plot waveforms stacked by depth
% figure; 
% hold on;
% 
% t = (1:nsamp) / sr * 1000;   % time in ms
% offset = max(abs(mean_wave_this_cluster(:))) * 1.5;  % vertical spacing
% 
% for i = 1:numel(channelsToPlot)
%     ch = channelsToPlot(i);
%     y_plot = mean_wave_this_cluster(ch,:) + (i-1)*offset;
% 
%     if ch == peakCh
%         plot(t, y_plot, 'r', 'LineWidth', 2.5); % peak channel in red
%     else
%         plot(t, y_plot, 'k', 'LineWidth', 1.5);
%     end
% end
% 
% % Add formatting
% xlabel('Time (ms)', 'FontSize', 12);
% ylabel('Depth (superficial → deep)', 'FontSize', 12);
% title(sprintf('Cluster %d (ID: %d) | Peak Ch %d @ (%.0f, %.0f) μm', ...
%        clusterID_index, good_clusters(clusterID_index), peakCh, peakCh_xpos, peakCh_ypos), ...
%        'FontSize', 13);
% 
% % Label with actual channel numbers and depths
% yticks((0:numel(channelsToPlot)-1) * offset);
% yticklabels(arrayfun(@(ch) sprintf('Ch %d (x=%.0f, y=%.0f μm)', ch, xpos(ch), ypos(ch)), ...
%             channelsToPlot, 'UniformOutput', false));
% 
% % Add text annotation
% text(0.02, 0.98, sprintf('Same column (x=%.0f μm)\nSpatial window: ±%.0f μm', ...
%      peakCh_xpos, spatialWindow), ...
%      'Units', 'normalized', 'VerticalAlignment', 'top', ...
%      'FontSize', 10, 'BackgroundColor', 'white');
% 
% grid on;
% box on;
% hold off;
% 
% % Print channel info to console
% fprintf('\n--- Cluster %d (ID: %d) ---\n', clusterID_index, good_clusters(clusterID_index));
% fprintf('Peak channel: %d at (x=%.0f, y=%.0f) μm\n', peakCh, peakCh_xpos, peakCh_ypos);
% fprintf('Plotting %d channels from same column (x=%.0f μm):\n', length(channelsToPlot), peakCh_xpos);
% for i = 1:length(channelsToPlot)
%     ch = channelsToPlot(i);
%     fprintf('  Ch %3d: x=%4.0f, y=%5.0f μm (Δy=%.0f μm)%s\n', ...
%             ch, xpos(ch), ypos(ch), abs(ypos(ch)-peakCh_ypos), ...
%             char(strcmp(num2str(ch), num2str(peakCh))*' <-- PEAK'));
% end
% fprintf('\n');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Visualize mean waveforms with spatial information preserved
% 
% clusterID_index = 98;  % which cluster INDEX (not cluster ID!)
% mean_wave_this_cluster = meanWav(:,:,clusterID_index);
% 
% % Find channel with largest spike amplitude
% [~, peakCh] = max(max(abs(mean_wave_this_cluster),[],2));
% 
% % Get the spatial position of the peak channel
% peakCh_ypos = ypos(peakCh);
% peakCh_xpos = xpos(peakCh);
% 
% % Get channels on the same shank
% sameShank = (xpos == peakCh_xpos);
% 
% % Calculate distance from peak channel (in y-direction)
% y_distance = abs(ypos - peakCh_ypos);
% 
% % Start with distance window
% spatialWindow = 150; % micrometers - adjust based on your needs
% nearby = sameShank & (y_distance <= spatialWindow);
% channelsToPlot = find(nearby);
% 
% % Sort by depth
% [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
% channelsToPlot = channelsToPlot(sortIdx);
% 
% % If too few channels, expand to nearest neighbors
% if length(channelsToPlot) < 11
%     warning('Only %d channels within %d μm, expanding to nearest neighbors', ...
%             length(channelsToPlot), spatialWindow);
% 
%     sameshank_channels = find(sameShank);
%     [~, sortIdx] = sort(y_distance(sameshank_channels));
%     sorted_channels = sameshank_channels(sortIdx);
% 
%     peakIdx = find(sorted_channels == peakCh);
%     chStart_idx = max(1, peakIdx - 5);
%     chEnd_idx = min(length(sorted_channels), peakIdx + 5);
%     channelsToPlot = sorted_channels(chStart_idx:chEnd_idx);
% 
%     [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
%     channelsToPlot = channelsToPlot(sortIdx);
% end
% 
% % If too many, keep closest to peak
% if length(channelsToPlot) > 11
%     [~, distIdx] = sort(y_distance(channelsToPlot));
%     channelsToPlot = channelsToPlot(distIdx(1:11));
% 
%     [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
%     channelsToPlot = channelsToPlot(sortIdx);
% end
% 
% % Plot waveforms stacked by depth
% figure; 
% % set(gcf, 'Position', [100, 100, 800, 600]);
% hold on;
% 
% t = (1:nsamp) / sr * 1000;   % time in ms
% offset = max(abs(mean_wave_this_cluster(:))) * 1.5;  % vertical spacing
% 
% for i = 1:numel(channelsToPlot)
%     ch = channelsToPlot(i);
%     y = mean_wave_this_cluster(ch,:) + (i-1)*offset;
% 
%     if ch == peakCh
%         plot(t, y, 'r', 'LineWidth', 2.5); % peak channel in red
%     else
%         plot(t, y, 'k', 'LineWidth', 1.5);
%     end
% end
% 
% % Add formatting
% xlabel('Time (ms)', 'FontSize', 12);
% ylabel('Depth (superficial → deep)', 'FontSize', 12);
% title(sprintf('Cluster %d (ID: %d) | Peak Ch %d @ %.0f μm', ...
%        clusterID_index, good_clusters(clusterID_index), peakCh, peakCh_ypos), ...
%        'FontSize', 13);
% 
% % Label with actual channel numbers and depths
% yticks((0:numel(channelsToPlot)-1) * offset);
% yticklabels(arrayfun(@(ch) sprintf('Ch %d (%.0f μm)', ch, ypos(ch)), ...
%             channelsToPlot, 'UniformOutput', false));
% 
% % Add vertical line at spike time
% plot([sbefore/sr*1000, sbefore/sr*1000], ylim, 'b--', 'LineWidth', 1);
% 
% % Add text annotation showing spatial window used
% if length(find(sameShank & (y_distance <= spatialWindow))) >= 11
%     text(0.02, 0.98, sprintf('Spatial window: ±%.0f μm', spatialWindow), ...
%          'Units', 'normalized', 'VerticalAlignment', 'top', ...
%          'FontSize', 10, 'BackgroundColor', 'white');
% else
%     text(0.02, 0.98, sprintf('Nearest neighbors (<%d channels in ±%.0f μm)', ...
%          length(find(sameShank & (y_distance <= spatialWindow))), spatialWindow), ...
%          'Units', 'normalized', 'VerticalAlignment', 'top', ...
%          'FontSize', 10, 'BackgroundColor', 'white', 'Color', 'red');
% end
% 
% grid on;
% box on;
% hold off;
% 
% % Optional: print channel info to console
% fprintf('\n--- Cluster %d (ID: %d) ---\n', clusterID_index, good_clusters(clusterID_index));
% fprintf('Peak channel: %d at %.0f μm\n', peakCh, peakCh_ypos);
% fprintf('Plotting %d channels:\n', length(channelsToPlot));
% for i = 1:length(channelsToPlot)
%     ch = channelsToPlot(i);
%     fprintf('  Ch %d: %.0f μm (%.0f μm from peak)%s\n', ...
%             ch, ypos(ch), abs(ypos(ch)-peakCh_ypos), ...
%             strcmp(ch, peakCh)*' <-- PEAK');
% end
% fprintf('\n');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% visualize mean waveforms (+/- 5 channels around the peak)
% 
% clusterID = 32;  % which cluster
% mean_wave_this_cluster = meanWav(:,:,clusterID);  % [nCh x nsamp]
% 
% % Find channel with largest spike
% [~, peakCh] = max(max(abs(mean_wave_this_cluster),[],2));
% 
% % Get +/- 5 channels around peak
% nChannelsAround = 5;
% chStart = max(1, peakCh - nChannelsAround);
% chEnd   = min(nCh, peakCh + nChannelsAround);
% 
% channelsToPlot = chStart:chEnd;
% 
% % Plot waveforms!
% figure; hold on;
% 
% t = (1:nsamp) / sr * 1000;   % time in ms
% offset = max(abs(mean_wave_this_cluster(:))) * 1.5;  % vertical spacing
% 
% for i = 1:numel(channelsToPlot)
%     ch = channelsToPlot(i);
%     y = mean_wave_this_cluster(ch,:) + (numel(channelsToPlot)-i)*offset;
% 
%     if ch == peakCh
%         plot(t, y, 'k', 'LineWidth', 2); % peak channel
%     else
%         plot(t, y, 'k', 'LineWidth', 1);
%     end
% end
% 
% xlabel('Time (ms)');
% ylabel('Channel (stacked)');
% title(['Cluster ' num2str(good_clusters(clusterID)) ...
%        ', channels ' num2str(chStart) '–' num2str(chEnd)]);
% 
% yticks((0:numel(channelsToPlot)-1) * offset);
% yticklabels(fliplr(channelsToPlot));
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % sort each cluster by probe depth
% clusterDepth = zeros(nClusters,1);
% clusterShank = zeros(nClusters,1);
% 
% for c = 1:nClusters
%     mw = meanWav(:,:,c);
%     [trough, ~] = min(mw,[],2);
%     amp = abs(trough);
% 
%     clusterDepth(c) = sum(amp .* ypos) / sum(amp);
% 
%     [~, maxCh] = max(amp);
%     clusterShank(c) = xpos(maxCh);
% end
% 
% % Shank-aware sorting
% uniqueShanks = unique(clusterShank);
% newOrder = [];
% for sh = uniqueShanks'
%     idx = find(clusterShank == sh);
%     [~, sidx] = sort(clusterDepth(idx),'descend'); % superficial → deep
%     newOrder = [newOrder; idx(sidx)];
% end
% 
% % Apply new order
% meanWav = meanWav(:,:,newOrder);
% good_clusters = good_clusters(newOrder);
% clusterDepth = clusterDepth(newOrder);
% clusterShank = clusterShank(newOrder);












