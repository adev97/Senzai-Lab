%% 1. Build geometry from the IMRO (field order: electrode, shank, bank, ref, channel)
imroFile = "R:\Basic_Sciences\Phys\SenzaiLab\Aparna\Shared_IMRO_Files\imro\forPaS_2Shank\shank3-4_staggered_dense.imro";   % full path if not in the current folder
txt = fileread(imroFile);
tok = regexp(txt, '\((\d+) (\d+) (\d+) (\d+) (\d+)\)', 'tokens');
M   = cell2mat(cellfun(@(c) str2double(c), tok(:), 'UniformOutput', false));
e = M(:,1);  shank = M(:,2);  ch = M(:,5) + 1;

xposIM = nan(384,1);  yposIM = nan(384,1);
xposIM(ch) = shank*250 + mod(e,2)*32;
yposIM(ch) = floor(e/2)*15;

% Kilosort positions (x should match the IMRO, y is the question)
chanPos = readNPY(fullfile(ksDir,'channel_positions.npy'));
xpos = chanPos(:,1);  ypos = chanPos(:,2);
fprintf('x identical: %d | y identical: %d\n', isequal(xpos,xposIM), isequal(ypos,yposIM));

%% 2. Load 20 s of raw data, filter, and common-median reference
nCh = 384;
fid = fopen(rawFile,'r');
fseek(fid, 60*sr*nCh*2, 'bof');               % start 60 s in; adjust if the file is shorter
raw = fread(fid, [nCh, sr*20], 'int16=>double');
fclose(fid);

raw = raw - mean(raw,2);
[b,a] = butter(3, 300/(sr/2), 'high');
f = filtfilt(b,a,raw')';
sd = std(f,0,2);
bad = sd > 5*median(sd) | sd < 0.2*median(sd);   % drop saturated/dead channels
f = f - median(f(~bad,:),1);
fprintf('%d channels flagged bad\n', sum(bad));

%% 3. Which geometry puts correlated channels next to each other?
rIM = adjCorr(f, xposIM, yposIM, bad);
rKS = adjCorr(f, xpos,   ypos,   bad);
fprintf('Adjacent-site correlation | IMRO: %.3f | Kilosort: %.3f\n', mean(rIM,'omitnan'), mean(rKS,'omitnan'));

%% 4. Where do the two geometries disagree? (24-row blocks, 360 um)
figure;
subplot(1,2,1); scatter(xposIM, yposIM, 12, 1:384, 'filled'); title('IMRO'); xlabel('x'); ylabel('y (um)'); colorbar;
subplot(1,2,2); scatter(xpos,   ypos,   12, 1:384, 'filled'); title('Kilosort'); xlabel('x'); colorbar;
sgtitle('Color = channel index. Smooth gradient = channels are physically contiguous');

%% local function
function r = adjCorr(f, x, y, bad)
r = [];
for s = unique(floor(x/250))'
    c = find(floor(x/250) == s & ~bad);
    [~,o] = sort(y(c));  c = c(o);
    for k = 1:numel(c)-1
        if y(c(k+1)) - y(c(k)) <= 15            % only truly adjacent rows
            r(end+1) = corr(f(c(k),:)', f(c(k+1),:)'); %#ok<AGROW>
        end
    end
end
end