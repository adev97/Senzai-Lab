% Example: load raw data

sr = 30000; % Probe sampling rate (Hz)
nchan_probe = 385;
chan = 1; % channel to load

% File path
% fp = "D:\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_FullFieldGrating\Record Node 102\experiment1\recording1\continuous\OneBox-100.ProbeA";

fp = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\SenzaiLab\Elissa_Belluccini\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_FullFieldGrating\Record Node 102\experiment1\recording1\continuous\OneBox-100.ProbeA";
cd(fp) 

d = LoadBinary('continuous.dat','frequency',sr,'nChannels',nchan_probe,...
'channels',chan);
d_dur = length(d)/sr; % Recording duration (sec)
td = (1:length(d))/sr; % Time array (sec)


% plot all raw data, entire timeseries (1600s)
figure; 
plot(td, d);
xlabel('Time (s)');
ylabel('Amplitude');
title('Raw Channel Data');   

%% filter raw data
% butterworth band pass filtering
low_cut = 300;      % Hz, removes LFP
high_cut = 6000;    % Hz, removes high-frequency noise
order = 3;          % 3rd-order filter (common for spike data)

% Design filter
[b, a] = butter(order, [low_cut high_cut] / (sr/2), 'bandpass');

% Apply filter (zero-phase to avoid shifting spikes)
d_filt = filtfilt(b, a, double(d));  % convert to double if not already

t0 = 10;          % start at 10 s
win = 0.05;       % 50 ms window
td_idx = td > t0 & td < (t0+win);

% plot filtered data on same original time axis
figure;
% plot(td, d_filt) %% plots entire timeseries
plot(td(td_idx), d_filt(td_idx));
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Filtered signal (300–6000 Hz)');


%% isolate spikes at timepoints?
% Robust noise estimate
noise_std = median(abs(d_filt)) / 0.6745;  % in µV

% Set threshold for spike detection
thr = -4 * noise_std;  % negative threshold for extracellular spikes

% detect threshold crossings
spike_idx = find(d_filt < thr);  % find points crossing negative threshold

% set refractory period
refractory = round(0.001 * sr);  % 1 ms in samples

% spike_times = spike_idx([true; diff(spike_idx) > refractory]);
% or
refractory = round(0.001 * sr);  % 1 ms
spike_times = [];

i = 1;
while i <= length(spike_idx)
    % find the contiguous points below threshold
    j = i;
    while j < length(spike_idx) && spike_idx(j+1) - spike_idx(j) == 1
        j = j + 1;
    end
    % take the minimum (negative peak) as spike time
    [~, min_idx] = min(d_filt(spike_idx(i:j)));
    spike_times(end+1,1) = spike_idx(i + min_idx - 1);

    % skip points within refractory period
    i = j + 1;
end

%% extract waveforms (1-2ms)
pre  = round(0.001 * sr);  % 1 ms before spike
post = round(0.002 * sr);  % 2 ms after spike

nSpikes = numel(spike_times);
waveforms = nan(nSpikes, pre + post + 1);

count = 0;
for i = 1:nSpikes
    t = spike_times(i);
    if t-pre > 0 && t+post <= length(d_filt)
        count = count + 1;
        waveforms(count,:) = d_filt(t-pre:t+post);
    end
end

waveforms = waveforms(1:count,:);


%% plot waveforms
figure;
plot(waveforms(:, :)', 'Color', [0 0 0 0.25]);

plot(waveforms(1:min(100,end), :)', 'Color', [0 0 0 0.25]);
xlabel('Samples');
ylabel('Amplitude (\muV)');
title('Single-channel spike waveforms');

%% plot waveforms on a time x axis
t_ms = (-pre:post) / sr * 1000;

figure;
plot(t_ms, waveforms(1:min(100,end), :)', 'Color', [0 0 0 0.25]);
hold on;
plot(t_ms, mean(waveforms,1), 'r', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title('Single-channel spike waveforms (time axis)');
