function [peakTimes, peakIdx, info] = detect(times, values, minDistanceMs)
%DETECT  spectHR's R-peak detector.
%   [PEAKTIMES, PEAKIDX, INFO] = Rpeak.detect(TIMES, VALUES, MINDISTANCEMS)
%   are the R-peaks of an upright, high-passed ECG, as
%   spectHR.signal.rpeak.detect_rpeaks finds them:
%     - the rate is one over the mean time step;
%     - a peak is a local maximum at least the median plus 1.5 standard
%       deviations of the signal;
%     - of peaks closer than MINDISTANCEMS (a refractory period; spectHR's
%       default is 300 ms, at most 200 beats a minute), the higher stays;
%       the distance in samples is truncated, as Python's int() does;
%     - each peak's time is moved by up to half a sample towards the higher
%       of its two neighbours: (next - previous) / rate / (2 * the larger of
%       its two steps).
%   PEAKTIMES are in TIMES' time base (seconds), PEAKIDX the peaks' samples
%   (1-based). INFO has the rate, the distance in samples and the threshold.
%
%   See also RPEAK.HIGHPASS, RPEAK.POLARITY, RPEAK.FINDPEAKS.
    times = double(reshape(times, 1, []));
    values = double(reshape(values, 1, []));
    peakTimes = zeros(1, 0);
    peakIdx = zeros(1, 0);
    info = struct('srate', NaN, 'minDistanceSamples', NaN, 'threshold', NaN);
    if numel(times) < 2 || numel(values) < 2
        return;
    end
    steps = diff(times);
    steps = steps(steps > 0);
    if isempty(steps)
        error('Rpeak:detect', 'The sampling rate cannot be told from these times.');
    end

    fs = 1 / Rpeak.npMean(steps);
    distance = max(1, fix((minDistanceMs / 1000) * fs));
    threshold = Rpeak.npMedian(values) + 1.5 * Rpeak.npStd(values);
    info = struct('srate', fs, 'minDistanceSamples', distance, 'threshold', threshold);

    peakIdx = Rpeak.findPeaks(values, 'Height', threshold, 'Distance', distance);
    if isempty(peakIdx)
        return;
    end
    n = numel(values);
    before = values(max(peakIdx - 1, 1));
    after = values(min(peakIdx + 1, n));
    at = values(peakIdx);
    contrast = max(abs(at - before), abs(after - at));
    contrast(contrast == 0) = 1e-12;
    peakTimes = times(peakIdx) + (after - before) / fs ./ (2 * contrast);
end
