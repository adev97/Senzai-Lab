function [trace_norm, valid, reason] = normalize_ccg(trace, flankIdx, minBaseline, minRefSpikes, nRef)
    minFlankSpikes = 5;
    valid = true; reason = '';

    if nRef < minRefSpikes
        valid = false;
        reason = sprintf('too few ref spikes (%d)', nRef);
        trace_norm = nan(size(trace));
        return
    end

    if sum(~isnan(trace(flankIdx))) < minFlankSpikes
        valid = false;
        reason = 'too few flank bins';
        trace_norm = nan(size(trace));
        return
    end

    baseline = mean(trace(flankIdx), 'omitnan');

    if baseline < minBaseline
        valid = false;
        reason = sprintf('baseline too low (%.2f spk/s)', baseline);
        trace_norm = nan(size(trace));
        return
    end

    trace_norm = trace / baseline;
