%% New Script that modifies DepthSort_meanWaveForms (from Yuta) for current Neuropixel data and only makes average waveform, without sorting

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("D:\buzcode-master")) 

sr = 30000;
nchan_probe = 384;
dtype = 'int16';

ksDir = 'D:\Kilosort\Mouse08_SC_20251007_810to2250\kilosort4';

spikeTimes    = readNPY(fullfile(ksDir,'spike_times.npy'));
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
nCh = size(chanPos,1);
dtype = 'int16';
bytesPerSample = 2;

fileInfo = dir(rawFile);
nSamples = fileInfo.bytes / (bytesPerSample * nCh);

sbefore = 30;   % samples before trough (-1 ms)
safter  = 60;   % samples after trough (+2 ms)
nsamp   = sbefore + safter;

maxSpikes = 1000;   % to sample per cluster

% pre-allocate
meanWav = zeros(nCh, nsamp, nClusters);

fid = fopen(rawFile,'r');

for c = 1:nClusters
    cluID = good_clusters(c);
    idx = find(spikeClusters == cluID);

    if numel(idx) > maxSpikes
        idx = idx(round(linspace(1,numel(idx),maxSpikes)));
    end

    wav = zeros(nCh, nsamp, numel(idx));

    for s = 1:numel(idx)
        t = spikeTimes(idx(s));

        if t <= sbefore || t + safter > nSamples
            continue
        end

        % Determine which chunk the spike is in
        fseek(fid, (t - sbefore - 1) * nCh * bytesPerSample, 'bof');
        w = fread(fid, [nCh, nsamp], '*int16');
        w = w - median(w(:,1:sbefore),2); % baseline
        wav(:,:,s) = w;
    end

    meanWav(:,:,c) = mean(wav,3);
end

fclose(fid);

% now you have the mean wave form for each good unit!

%% visualize mean waveforms (+/- 5 channels around the peak)












% sort each cluster by probe depth
clusterDepth = zeros(nClusters,1);
clusterShank = zeros(nClusters,1);

for c = 1:nClusters
    mw = meanWav(:,:,c);
    [trough, ~] = min(mw,[],2);
    amp = abs(trough);

    clusterDepth(c) = sum(amp .* ypos) / sum(amp);

    [~, maxCh] = max(amp);
    clusterShank(c) = xpos(maxCh);
end

% Shank-aware sorting
uniqueShanks = unique(clusterShank);
newOrder = [];
for sh = uniqueShanks'
    idx = find(clusterShank == sh);
    [~, sidx] = sort(clusterDepth(idx),'descend'); % superficial → deep
    newOrder = [newOrder; idx(sidx)];
end

% Apply new order
meanWav = meanWav(:,:,newOrder);
good_clusters = good_clusters(newOrder);
clusterDepth = clusterDepth(newOrder);
clusterShank = clusterShank(newOrder);























% load raw data file (Merged.dat)
rawFile = "D:\Kilosort\Mouse08_SC_20251007_810to2250\Merged.dat";
nCh = size(chanPos,1);
dtype = 'int16';

fileInfo = dir(rawFile);
bytesPerSample = 2;
nSamples = fileInfo.bytes / (bytesPerSample * nCh);

mm = memmapfile(rawFile, 'Format', {'int16', [nCh, nSamples], 'data'});

sbefore = 20;   % samples before trough (~0.67 ms)
safter  = 44;   % samples after trough (~1.47 ms)
nsamp   = sbefore + safter;

maxSpikes = 1000;   % to sample per cluster

% mean waveform per cluster
nClusters = numel(good_clusters);
meanWav = zeros(nCh, nsamp, nClusters);
clusterDepth = zeros(nClusters,1);
clusterShank = zeros(nClusters,1);

for c = 1:nClusters
    cluID = good_clusters(c);
    idx = find(spikeClusters == cluID);

    if numel(idx) > maxSpikes
        idx = idx(round(linspace(1,numel(idx),maxSpikes)));
    end

    wav = zeros(nCh, nsamp, numel(idx));

    for s = 1:numel(idx)
        if t <= sbefore || t + safter > nSamples
            continue
        end

        t = spikeTimes(idx(s));
        w = mm.data(:, t - sbefore + 1 : t + safter);
        w = w - median(w(:,1:sbefore),2);   % baseline
        wav(:,:,s) = w;
    end

    mw = mean(wav,3);
    meanWav(:,:,c) = mw;

    [trough, ~] = min(mw,[],2);
    amp = abs(trough);

    clusterDepth(c) = sum(amp .* ypos) / sum(amp);

    [~, maxCh] = max(amp);
    clusterShank(c) = xpos(maxCh);
end


% depth sorting based on which shank it is
uniqueShanks = unique(clusterShank);

newOrder = [];
for sh = uniqueShanks'
    idx = find(clusterShank == sh);
    [~, sidx] = sort(clusterDepth(idx),'descend'); % superficial → deep
    newOrder = [newOrder; idx(sidx)];
end













waveforms = zeros(length(good_clusters), nchan_probe, 1000); % Preallocate for waveforms
for i = 1:length(good_clusters)
    clusterSpikes = spikeTimes(spikeClusters == good_clusters(i));
    for j = 1:length(clusterSpikes)
        % Extract waveforms around each spike time
        waveforms(i, :, j) = templates(:, good_clusters(i), 1); % Assuming templates are organized by cluster
    end
end

meanWaveforms = squeeze(mean(waveforms, 3)); % Average waveforms across spikes


