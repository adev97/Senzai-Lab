%% Script to print which probe channel (ON MAP FROM PYTHON) is associated with which recording index

% Load the channel map
chanMap = readNPY(fullfile(ksDir,'channel_map.npy'));

% Check what probe channels correspond to recording indices
fprintf('Recording index 341 → Probe channel %d\n', chanMap(341));
fprintf('Recording index 287 → Probe channel %d\n', chanMap(287));

% Find which recording indices are on the same shank as recording index 341
peakCh_recording_idx = 341;
peakCh_xpos = xpos(peakCh_recording_idx);
sameshank_recording_idx = find(xpos == peakCh_xpos);

fprintf('\nRecording indices on same shank as 341:\n');
fprintf('%d ', sameshank_recording_idx);
fprintf('\n');

fprintf('\nTheir corresponding probe channels:\n');
for i = 1:min(10, length(sameshank_recording_idx))
    fprintf('Recording idx %d → Probe ch %d\n', ...
        sameshank_recording_idx(i), chanMap(sameshank_recording_idx(i)));
end