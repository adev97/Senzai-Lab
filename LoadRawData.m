% Example: load raw data

sr = 30000; % Probe sampling rate (Hz)
nchan_probe = 385;
chan = 1; % channel to load

% File path
fp = "D:\OpenEphys_Data\Mouse08\Mouse08_20251007_4shanks_810to2250_FullFieldGrating\Record Node 102\experiment1\recording1\continuous\OneBox-100.ProbeA";
cd(fp) 

d = LoadBinary('continuous.dat','frequency',sr,'nChannels',nchan_probe,...
'channels',chan);
d_dur = length(d)/sr; % Recording duration (sec)
td = (1:length(d))/sr; % Time array (sec)
   