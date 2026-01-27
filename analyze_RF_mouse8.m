%% Analyze unit responses to specific stimuli (receptive field)

addpath(genpath("C:\Users\urs2027\Documents\GitHub\Senzai-Lab\buzcode-master"));

sr = 30000; % Hz
nchan_probe = 384;
dtype = 'int16';

%% get spike times per stimulus (from events folder) 
% use TTLs to align spike onset
TTL_OneBox = "D:\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_RFMapping\Record Node 102\experiment1\recording1\events\OneBox-100.OneBox-ADC\TTL"; 

eventSamples = readNPY(fullfile(TTL_OneBox, 'sample_numbers.npy')); 
states = readNPY(fullfile(TTL_OneBox,'states.npy')); 
fullWords = readNPY(fullfile(TTL_OneBox,'full_words.npy'));

stimOn_idx = states == 1;
stimOn_samples = eventSamples(stimOn_idx);
stimWords = fullWords(stimOn_idx);

pre_ms = 50;
post_ms = 200;

pre_samp  = round(pre_ms * sr / 1000);
post_samp = round(post_ms * sr / 1000);

unitID = good_clusters(1);
spk = spikeTimes(spikeClusters == unitID);

allSpikesAligned = [];

for k = 1:length(stimOn_samples)
    t0 = stimOn_samples(k);
    rel_spk = spk(spk >= t0-pre_samp & spk <= t0+post_samp) - t0;
    allSpikesAligned = [allSpikesAligned; rel_spk];
end

histogram(allSpikesAligned / sr * 1000, -50:2:200);
xlabel('Time from stim onset (ms)');
ylabel('Spike count');
title(sprintf('Unit %d PSTH (all stimuli)', unitID));
