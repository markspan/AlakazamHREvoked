function [beatTimes, ibiMs, polarityInfo, nStretches] = beats(values, timing, breakAfter, opt)
%BEATS  spectHR's whole route on one ECG signal: the R-peaks, their times,
%   and the interval from each to the next.
%
%   [BEATTIMES, IBIMS, POLARITYINFO, NSTRETCHES] = Rpeak.beats(VALUES, TIMING,
%   BREAKAFTER, OPT) finds the R-peaks of the ECG VALUES as spectHR does:
%     1. its polarity once, on all its valid samples (Rpeak.polarity), or as
%        OPT.Polarity says ('Upright', 'Inverted'); an inverted ECG is turned
%        upright to search it;
%     2. per stretch, the high-pass (Rpeak.highpass, OPT.HighpassHz, 0 for
%        none) and the detector (Rpeak.detect, OPT.MinIntervalMs).
%
%   TIMING is the sampling rate of a regular signal (an EEG channel), whose
%   stretches are then timed from their own start, (0:m-1)/rate, as spectHR
%   times a recording of its own, and shifted after; or the time of every
%   sample (a stream's timestamps, seconds), used as they are.
%
%   A STRETCH is a run of valid (not NaN) samples that BREAKAFTER (a logical
%   per sample, or []) does not cut, nor a gap in the times longer than 2.5
%   sample steps. No beat is found across one, and no interval spans one.
%
%   BEATTIMES are in TIMING's time base: seconds from the first sample for a
%   rate, the stream's own clock for times. IBIMS is each beat's interval to
%   the next in ms, NaN for the last of a stretch. POLARITYINFO has polarity
%   ('normal' or 'inverted'), decidedBy and skewness.
%
%   See also RPEAK.DETECT, RPEAK.HIGHPASS, RPEAK.POLARITY, HEARTBEATS.
    values = double(reshape(values, 1, []));
    n = numel(values);
    regular = isscalar(timing);
    if regular
        srate = double(timing);
        times = (0:n - 1) / srate;
    else
        times = double(reshape(timing, 1, []));
    end
    if isempty(breakAfter)
        breakAfter = false(1, n);
    end
    breakAfter = reshape(logical(breakAfter), 1, []);
    if ~regular && n > 1
        steps = diff(times);
        breakAfter(1:n - 1) = breakAfter(1:n - 1) | steps > 2.5 * median(steps(steps > 0));
    end

    valid = ~isnan(values);
    polarityInfo = struct('polarity', 'normal', 'decidedBy', 'the user', 'skewness', NaN);
    switch lower(char(string(opt.Polarity)))
        case 'upright'
        case 'inverted'
            polarityInfo.polarity = 'inverted';
        otherwise
            % A regular signal's valid samples timed as one run, as before
            % 0.2.0; a stream's by their own times.
            if regular
                polarityTimes = (0:nnz(valid) - 1) / srate;
            else
                polarityTimes = times(valid);
            end
            [polarityInfo.polarity, info] = Rpeak.polarity(polarityTimes, values(valid));
            polarityInfo.skewness = info.skewness;
            polarityInfo.decidedBy = strrep(info.decidedBy, '_', ' ');
    end
    if strcmp(polarityInfo.polarity, 'inverted')
        values = -values;
    end

    starts = find(valid & [true, ~valid(1:end - 1) | breakAfter(1:end - 1)]);
    ends = find(valid & [~valid(2:end) | breakAfter(1:end - 1), true]);
    nStretches = numel(starts);
    beatTimes = zeros(1, 0);
    ibiMs = zeros(1, 0);
    for s = 1:nStretches
        span = starts(s):ends(s);
        if regular
            t = (0:numel(span) - 1) / srate;
        else
            t = times(span);
        end
        x = values(span);
        if opt.HighpassHz > 0
            x = Rpeak.highpass(x, t, opt.HighpassHz);
        end
        found = Rpeak.detect(t, x, opt.MinIntervalMs);
        if regular
            found = found + (starts(s) - 1) / srate;
        end
        beatTimes = [beatTimes, found]; %#ok<AGROW>
        ibiMs = [ibiMs, diff(found) * 1000, NaN(1, min(1, numel(found)))]; %#ok<AGROW>
    end
end
