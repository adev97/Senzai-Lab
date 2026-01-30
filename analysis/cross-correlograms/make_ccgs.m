%% Make cross-correlograms (CCGs) for different units to assess likelyhood of firing together

% CCG.m from buzcode

addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master")) 

% example function usage from yuta
% [ccg,t_ccg] = CCG(ts_res,gs_res,'binSize',0.02,'duration',8);

% goal: identify units that are firing together

sr = 30000;
nchan_probe = 384;
dtype = 'int16';

% kilosort directory
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

%% spikeTimes - contains all spike timings (when)
%% spikeClusters - contains which unit gave the spike (who) and preserves channel number structure

%% Now building cross-correlograms, searching for units that fire together
% during this data collection

% convert time to seconds instead of samples
ts = double(spikeTimes) / sr;   % convert to seconds
gs = double(spikeClusters);     % unit IDs

% Map original cluster IDs to consecutive indices
[~, gs_idx] = ismember(spikeClusters, good_clusters);

% Use gs_idx for CCG
[ccg, t] = CCG(ts, gs_idx, 'binSize', binSize, 'duration', duration);


