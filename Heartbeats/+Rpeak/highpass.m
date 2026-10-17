function y = highpass(values, times, cutoff)
%HIGHPASS  spectHR's ECG prefilter: a zero-phase order-4 Butterworth high-pass.
%   Y = Rpeak.highpass(VALUES, TIMES, CUTOFF) is what spectHR's
%   Samples.filtered(filter_type='highpass', cutoff=CUTOFF) gives (order 4,
%   its default there), the step before beat detection that strips baseline
%   wander (spectHR's default cutoff is 0.5 Hz). Its rate is spectHR's
%   Samples.srate, one over the median time step. Like spectHR, a signal
%   too short for the filter is returned unfiltered.
%
%   See also RPEAK.BUTTER, RPEAK.FILTFILT, RPEAK.DETECT.
    values = double(reshape(values, 1, []));
    srate = 1 / Rpeak.npMedian(diff(double(times)));
    [b, a] = Rpeak.butter(4, cutoff / (0.5 * srate), 'high');
    if numel(values) <= 3 * max(numel(a), numel(b))
        y = values;     % spectHR's filter_ecg: the prefilter skipped
        return;
    end
    y = Rpeak.filtfilt(b, a, values);
end
