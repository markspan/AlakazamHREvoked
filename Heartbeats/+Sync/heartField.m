function check = heartField(EEG, latency, exclude)
%HEARTFIELD  The beats' timing checked on the heart's own field in the EEG.
%
%   CHECK = Sync.heartField(EEG, LATENCY, EXCLUDE) averages the scalp EEG at
%   the beats' latencies (samples) and finds where the heart's electrical
%   field peaks around them. The heart's field reaches every scalp electrode,
%   largest during the QRS complex (Dirlich et al., 1997); it is there at
%   the same moment as the R-wave in the ECG, so the peak should fall at
%   0 ms, and where it falls instead is the timing's error from end to end,
%   whatever its cause: clocks not lined up, or a device that delivers its
%   samples late (a delay within a device, which Lab Streaming Layer cannot
%   measure).
%
%   The field's strength over time is the global field power (the standard
%   deviation across the scalp channels) of the beat-locked average, from
%   200 ms before each beat to 200 ms after, less each channel's mean from
%   -200 to -120 ms. Its peak within 100 ms of the beats, placed to a
%   fraction of a sample (a parabola through the peak and its neighbours),
%   is the lag. The field counts as seen when that peak is at least three
%   times the median power before the beats; otherwise nothing is concluded.
%   The QRS of the heart's field and the R-wave of one ECG lead need not peak
%   at quite the same moment, so the check confirms the timing to within a
%   few ms: it catches what is wrong by more.
%
%   EXCLUDE are channels to leave out (the ECG channel itself). CHECK has
%   seen, lagMs, ratio (the peak over the power before the beats), beats
%   and channels (how many were averaged).
%
%   See also HEARTBEATS.
    check = struct('seen', false, 'lagMs', NaN, 'ratio', NaN, 'beats', 0, 'channels', 0);
    srate = double(EEG.srate);
    channels = find(reshape(logical(eegChannelMask(EEG.chanlocs)), 1, []));
    channels = setdiff(channels, exclude);
    check.channels = numel(channels);
    if numel(channels) < 2 || isempty(latency)
        return;
    end

    offsets = -round(0.2 * srate):round(0.2 * srate);
    centres = round(latency);
    centres = centres(centres + offsets(1) >= 1 & centres + offsets(end) <= size(EEG.data, 2));
    data = double(EEG.data(channels, :));
    total = zeros(numel(channels), numel(offsets));
    for c = centres
        window = data(:, c + offsets);
        if ~any(isnan(window(:)))
            total = total + window;
            check.beats = check.beats + 1;
        end
    end
    if check.beats < 20
        return;
    end

    average = total / check.beats;
    t = offsets / srate * 1000;
    before = t <= -120;
    average = average - mean(average(:, before), 2);
    power = std(average, 1, 1);
    search = find(abs(t) <= 100);
    [peak, k] = max(power(search));
    at = search(k);
    shift = 0;
    if at > 1 && at < numel(power)
        curve = power(at - 1) - 2 * power(at) + power(at + 1);
        if curve < 0
            shift = 0.5 * (power(at - 1) - power(at + 1)) / curve;
        end
    end
    check.lagMs = t(at) + shift * 1000 / srate;
    check.ratio = peak / median(power(before));
    check.seen = check.ratio >= 3;
end
