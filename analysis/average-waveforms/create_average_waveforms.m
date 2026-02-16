%% New Script that modifies DepthSort_meanWaveForms (from Yuta) for current Neuropixel data and only makes average waveform, without sorting

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master")) 

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
% cgFile = fullfile(ksDir,'cluster_group.tsv'); % pre-manual curation
cgFile = fullfile(ksDir, 'cluster_KSLabel.tsv'); % post-manual curation

cluster_groups = readtable(cgFile, ...
    'FileType','text', ...
    'Delimiter','\t');

keepGroups = {'good'}; % only keep good labeled units %%%%%%%%%%%%%%

% if using pre-manual cluster_group.tsv
% toKeep = ismember(cluster_groups.group, keepGroups);
% keepClusters = cluster_groups.cluster_id(toKeep); % contains only the good cluster ids (162 total)

% if using post-manual cluster_KSLabel.tsv
toKeep = ismember(cluster_groups.KSLabel, keepGroups);
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

% now you have the mean wave form for each good unit for each channel!

%% Visualize mean waveforms - filter by channel parity to avoid checkerboard (even vs odd units)

% changed folder structure, change outputs if need to run again
outWaveDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\20251007_MergedDat\cluster_KSLabel\plus2ms-minus1ms\average-waveforms-per-unit-from-raw';
outDepthDir = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\20251007_MergedDat\cluster_KSLabel\plus2ms-minus1ms\heatmap-21chan-per-unit-from-raw';

if ~exist(outWaveDir, 'dir'); mkdir(outWaveDir); end
if ~exist(outDepthDir, 'dir'); mkdir(outDepthDir); end


for clusterID_index = 55:nClusters

    mean_wave_this_cluster = meanWav(:,:,clusterID_index);
    num_channels_to_plot = 21; %% change this based on how many channels you want to visualize above/below the peak channel

    % Find channel with largest spike amplitude
    [~, peakCh] = max(max(abs(mean_wave_this_cluster),[],2));
    
    peakCh_ypos = ypos(peakCh);
    peakCh_xpos = xpos(peakCh);

    % Get channels on same shank (same x position)
    sameshank_channels = find(xpos == peakCh_xpos);

    % Filter by same parity (odd/even recording index) as peak channel
    peakParity = mod(peakCh, 2);  % 0 for even, 1 for odd
    sameParity = (mod(sameshank_channels, 2) == peakParity);
    sameParity_channels = sameshank_channels(sameParity);

    % Calculate distance from peak channel
    y_distance = abs(ypos(sameParity_channels) - peakCh_ypos);

    % Start with spatial window
    spatialWindow = 300; % micrometers
    nearby_idx = find(y_distance <= spatialWindow);
    channelsToPlot = sameParity_channels(nearby_idx);

    % If too few channels within window, expand to nearest neighbors
    if length(channelsToPlot) < num_channels_to_plot
        [~, sortIdx] = sort(y_distance);
        n_to_take = min(num_channels_to_plot, length(sameParity_channels));
        channelsToPlot = sameParity_channels(sortIdx(1:n_to_take));
    end

    % If too many channels, keep closest to peak
    if length(channelsToPlot) > num_channels_to_plot
        distances_to_plot = abs(ypos(channelsToPlot) - peakCh_ypos);
        [~, distIdx] = sort(distances_to_plot);
        channelsToPlot = channelsToPlot(distIdx(1:num_channels_to_plot));
    end

    % Sort by depth for plotting (superficial to deep)
    [~, sortIdx] = sort(ypos(channelsToPlot), 'descend');
    channelsToPlot = channelsToPlot(sortIdx);

    %% Waveforms stacked by depth
    fig1 = figure('Visible','off');
    hold on;
    
    t = (1:nsamp) / sr * 1000;   % time in ms
    offset = max(abs(mean_wave_this_cluster(:))) * 1.5;  % vertical spacing
    
    for i = 1:numel(channelsToPlot)
        ch = channelsToPlot(i);
        y_plot = mean_wave_this_cluster(ch,:) + (i-1)*offset;
        
        if ch == peakCh
            plot(t, y_plot, 'r', 'LineWidth', 2);
        else
            plot(t, y_plot, 'k', 'LineWidth', 1);
        end
    end
    
    xlabel('Time (ms)');
    ylabel('Depth (superficial → deep)');
    title(sprintf('Unit %d | Peak Ch %d @ %.0f μm', ...
            good_clusters(clusterID_index), peakCh, peakCh_ypos));
     
    peakIdxInPlot = find(channelsToPlot == peakCh);
    peakZeroY = (peakIdxInPlot - 1) * offset;
    yline(peakZeroY, ':', 'LineWidth', 0.75, 'Color', [0.4 0.4 0.4]);
    
    yticks((0:numel(channelsToPlot)-1) * offset);
    yticklabels(arrayfun(@(ch) sprintf('Ch %d (%.0f μm)', ch, ypos(ch)), ...
                channelsToPlot, 'UniformOutput', false));
    
    % SAVE WAVEFORMS
    saveas(fig1, fullfile(outWaveDir, ...
        sprintf('unit_%03d_avg.png', good_clusters(clusterID_index))));
    close(fig1);  

    
    %% Heatmap spike propogation
    
    fig2 = figure('Visible','off');
    % subplot(1,2,1);
    imagesc(t, 1:length(channelsToPlot), mean_wave_this_cluster(channelsToPlot,:));
    colormap('jet');
    colorbar;
    xlabel('Time (ms)');
    ylabel('Channel (superficial → deep)');
    title('Waveform heatmap');
    set(gca, 'YTick', 1:length(channelsToPlot));
    set(gca, 'YTickLabel', arrayfun(@(ch) sprintf('Ch%d', ch), channelsToPlot, 'UniformOutput', false));
    
    % Find peak time on each channel
    % subplot(1,2,2);
    % peak_times = zeros(length(channelsToPlot), 1);
    % for i = 1:length(channelsToPlot)
    %     ch = channelsToPlot(i);
    %     [~, peak_idx] = min(mean_wave_this_cluster(ch,:)); % trough time
    %     peak_times(i) = t(peak_idx);
    % end
    % 
    % plot(peak_times, ypos(channelsToPlot), 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    % xlabel('Trough time (ms)');
    % ylabel('Depth (μm)');
    % title('Spike propagation');
    % grid on;
    
    saveas(fig2, fullfile(outDepthDir, ...
        sprintf('unit_%03d_heatmap_21chan.png', good_clusters(clusterID_index))));
    close(fig2);
    
end

%% Save mean waveforms and cluster IDs for cell type identification (python script)
saveFolder = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Aparna\Mouse08\mean-waveforms-good-clusters';
save(fullfile(saveFolder, 'meanWav_units.mat'), 'meanWav', 'good_clusters', 'xpos', 'ypos');
disp("Mean Waveforms Saved!")








































