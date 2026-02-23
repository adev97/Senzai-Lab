function zJitterFull = computeJitterZ(s_state, candPairs, pairs, postWin, binSize, nShuffle, jitterWindow)

nCand = length(candPairs);
zJitter = nan(nCand,1);

jitter_spikes = @(ts) ts + (rand(size(ts)) - 0.5) * jitterWindow;

fprintf('Applying jitter correction on %d candidate pairs...\n', nCand);

parfor idx = 1:nCand
    p = candPairs(idx);
    i = pairs(p,1);
    j = pairs(p,2);

    % real peak
    ts_i = s_state(s_state(:,2)==i,1);
    ts_j = s_state(s_state(:,2)==j,1);

    s_tmp_real = [ts_i ones(size(ts_i)); ts_j 2*ones(size(ts_j))];
    [ccg_real, ~] = CCG(s_tmp_real(:,1), s_tmp_real(:,2), ...
                        'binSize', binSize, 'duration', 0.04);
    realPeak = sum(ccg_real(postWin,1,2));

    % shuffle peaks
    shufflePeaks = zeros(nShuffle,1);
    for sIter = 1:nShuffle
        ts_i_jit = jitter_spikes(ts_i);
        ts_j_jit = jitter_spikes(ts_j);

        s_tmp = [ts_i_jit ones(size(ts_i_jit));
                 ts_j_jit 2*ones(size(ts_j_jit))];

        [ccg_tmp, ~] = CCG(s_tmp(:,1), s_tmp(:,2), ...
                            'binSize', binSize, 'duration', 0.04);
        shufflePeaks(sIter) = sum(ccg_tmp(postWin,1,2));
    end

    % compute z-score
    muJ = mean(shufflePeaks);
    sdJ = std(shufflePeaks);
    if sdJ > 0
        zJitter(idx) = (realPeak - muJ)/sdJ;
    end
end

% assign back to full array
zJitterFull = nan(max(pairs(:)),1);
zJitterFull(candPairs) = zJitter;

fprintf('Jitter correction done.\n');

