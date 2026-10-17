function [polarity, info] = polarity(times, values)
%POLARITY  spectHR's ECG polarity check: is the ECG upright or inverted?
%   [POLARITY, INFO] = Rpeak.polarity(TIMES, VALUES) is 'normal' or
%   'inverted', as spectHR.signal.ecg.detect_ecg_polarity decides it for a
%   whole recording:
%     - the middle third of the signal (both ends settle), less its median;
%     - band-passed 5 to 20 Hz (order 3, zero phase), where the QRS
%       complexes dominate;
%     - the sign of its skewness: tall upward R-waves make it positive;
%     - only where the skewness is too small to trust (below 0.1 either
%       way), the summed prominence of the upward peaks against that of the
%       downward ones (peaks at least 250 ms apart and a standard deviation
%       prominent); a flat signal is taken as upright.
%   INFO has the skewness, the peak score, and what decided ('skewness',
%   'peak_count', or 'raw_skewness' for a signal too short to filter).
%
%   See also RPEAK.SKEW, RPEAK.FINDPEAKS, RPEAK.DETECT.
    times = double(reshape(times, 1, []));
    values = double(reshape(values, 1, []));
    steps = diff(times);
    fs = 1 / Rpeak.npMean(steps(steps > 0));

    n = numel(values);
    trim = floor(n / 3);
    ecg = values;
    if n > 2 * trim && trim > 0
        ecg = values(trim + 1:n - trim);
    end
    ecg = ecg - Rpeak.npMedian(ecg(~isnan(ecg)));

    info = struct('skewness', NaN, 'peakScore', 0, 'decidedBy', 'skewness');
    if numel(ecg) >= 50
        nyq = 0.5 * fs;
        band = [max(5 / nyq, 1e-4), min(20 / nyq, 1 - 1e-4)];
        [b, a] = Rpeak.butter(3, band, 'band');
        filtered = Rpeak.filtfilt(b, a, ecg);
    else
        filtered = ecg;
        info.decidedBy = 'raw_skewness';
    end

    info.skewness = Rpeak.skew(filtered);
    if abs(info.skewness) >= 0.1
        if info.skewness > 0
            polarity = 'normal';
        else
            polarity = 'inverted';
        end
        return;
    end

    info.decidedBy = 'peak_count';
    distance = max(fix(0.25 * fs), 1);
    floorProminence = Rpeak.npStd(filtered);
    [~, up] = Rpeak.findPeaks(filtered, 'Distance', distance, 'Prominence', floorProminence);
    [~, down] = Rpeak.findPeaks(-filtered, 'Distance', distance, 'Prominence', floorProminence);
    upward = 0;
    downward = 0;
    if ~isempty(up)
        upward = Rpeak.npSum(up);
    end
    if ~isempty(down)
        downward = Rpeak.npSum(down);
    end
    info.peakScore = upward - downward;
    if upward == 0 && downward == 0
        polarity = 'normal';
    elseif upward >= downward
        polarity = 'normal';
    else
        polarity = 'inverted';
    end
end
