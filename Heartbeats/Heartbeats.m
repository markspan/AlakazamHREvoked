function [EEG, options] = Heartbeats(input, varargin)
%% HEARTBEATS  Every heartbeat as an event, from the R-peaks of an ECG channel:
%   the first step to heartbeat-evoked potentials.
%
%   The R-peaks are found exactly as spectHR finds them, ported to MATLAB
%   (the +Rpeak package, held to spectHR's own results by
%   SpectHRParityTest):
%     1. the ECG's polarity, from the skewness of its 5 to 20 Hz band (an
%        inverted ECG is turned upright for detection, not in the data);
%     2. a zero-phase 0.5 Hz high-pass, order 4, against baseline wander;
%     3. the local maxima at least the median plus 1.5 standard deviations,
%        no two closer than the shortest interval (300 ms), each placed to a
%        fraction of a sample from its neighbours.
%   The data are not changed. Each R-peak becomes an event of the code
%   chosen ('heartbeat' by default), at its sub-sample latency, with two
%   fields: IBI, the interval to the next beat in ms (as spectHR's IBI
%   series), and HR, the heart rate it makes in beats per minute; NaN for
%   the last beat before a cut, a rejected stretch or the end.
%
%   THE BINS come next, from DefineBins, which can choose heartbeats by what
%   happens around them, e.g. the heartbeats in the second after a face:
%       bin 1 "HEP after faces" : "heartbeat" prev(1) within [-1000, 0] ms
%   and Average then gives the heartbeat-evoked potential of each bin.
%
%   CUTS AND REJECTED STRETCHES. A recording cut by 'boundary' events, or
%   with rejected (NaN) samples, is a set of separate stretches: each is
%   filtered and searched on its own, as spectHR would a recording of its
%   own, so no beat is found across a join. The polarity is decided once,
%   on all of the channel's valid samples.
%
%   THE ECG's ELECTRICAL FIELD in the EEG (the cardiac field artefact) is
%   not removed here: remove the heart components with ICA first.
%
%   WHAT IT RECORDS, in EEG.etc.alz.heartbeats: the plugin version, the
%   channel, the polarity and what decided it, the settings, the number of
%   beats and the mean heart rate.
%
%   Signature (Alakazam transformation contract):
%     [EEG, options] = Heartbeats(input)        % interactive dialog
%     [EEG, options] = Heartbeats(input, opts)  % replay a stored options struct
%
%   See also RPEAK.DETECT, RPEAK.POLARITY, RPEAK.HIGHPASS, DEFINEBINS.
[opts, interactive] = TransTools.InitGuard(nargin, 'Alakazam:Heartbeats', varargin{:});

requireContinuous(input);

if interactive
    stored = TransformSettings.get('Heartbeats');
    if isempty(stored) || ~isstruct(stored)
        stored = struct();
    end
    d = @(f, v) TransTools.FieldOr(stored, f, v);
    options = TransformOptionsDialog( ...
        'title', 'Heartbeats', ...
        'Description', ['Every R-peak of the ECG channel becomes an event, found as spectHR ' ...
            'finds it. The data are not changed. Choose the heartbeats per condition with ' ...
            'DefineBins afterwards, and Average them.'], ...
        {'ECG channel'; 'Channel'}, DialogFields.Channels(input, d('Channel', ecgChannelGuess(input)), ...
            'Multiple', false, 'Required', true, 'QuickSelect', cell(0, 2)), ...
        {'ECG polarity'; 'Polarity'}, DialogFields.Choice({'Automatic', 'Upright', 'Inverted'}, ...
            d('Polarity', 'Automatic')), ...
        {'Shortest interval (ms)'; 'MinIntervalMs'}, DialogFields.Number(d('MinIntervalMs', 300), ...
            'Limits', [1 Inf]), ...
        {'High-pass before detection (Hz)'; 'HighpassHz'}, DialogFields.Number(d('HighpassHz', 0.5), ...
            'Limits', [0 Inf]), ...
        {'Event code'; 'EventType'}, DialogFields.Text(d('EventType', 'heartbeat')));
    if isempty(options)
        EEG = [];   % cancelled: no node, no compute
        return;
    end
    TransformSettings.set('Heartbeats', options);
else
    options = opts;
end
opt = normaliseOptions(options, input);

EEG = input;
channel = find(strcmpi({EEG.chanlocs.labels}, opt.Channel), 1);
ecg = double(EEG.data(channel, :));
srate = double(EEG.srate);

[polarity, decidedBy, skewness] = choosePolarity(ecg, srate, opt.Polarity);
upright = ecg;
if strcmp(polarity, 'inverted')
    upright = -ecg;
end

latency = zeros(1, 0);
ibi = zeros(1, 0);
stretches = validStretches(EEG, ecg);
for s = 1:size(stretches, 1)
    first = stretches(s, 1);
    span = first:stretches(s, 2);
    times = (0:numel(span) - 1) / srate;
    values = upright(span);
    if opt.HighpassHz > 0
        values = Rpeak.highpass(values, times, opt.HighpassHz);
    end
    peakTimes = Rpeak.detect(times, values, opt.MinIntervalMs);
    latency = [latency, first + peakTimes * srate]; %#ok<AGROW>
    ibi = [ibi, diff(peakTimes) * 1000, NaN(1, min(1, numel(peakTimes)))]; %#ok<AGROW>
end
if isempty(latency)
    throw(MException('Alakazam:Heartbeats', ['I''m afraid no heartbeat was found in %s. ' ...
        'Is it the ECG channel, and is the ECG in it?'], opt.Channel));
end

[EEG, replaced] = addBeatEvents(EEG, latency, ibi, opt.EventType);
hr = 60000 ./ ibi;
EEG.etc.alz.heartbeats = struct( ...
    'version',        pluginVersion(), ...
    'channel',        opt.Channel, ...
    'polarity',       polarity, ...
    'polarityBy',     decidedBy, ...
    'skewness',       skewness, ...
    'minIntervalMs',  opt.MinIntervalMs, ...
    'highpassHz',     opt.HighpassHz, ...
    'eventType',      opt.EventType, ...
    'beats',          numel(latency), ...
    'stretches',      size(stretches, 1), ...
    'meanHr',         mean(hr(isfinite(hr))), ...
    'replaced',       replaced, ...
    'options',        options);

fprintf('Heartbeats: %d heartbeats in %s (%s ECG, by %s), %.1f beats a minute on average%s.\n', ...
    numel(latency), opt.Channel, upright2word(polarity), decidedBy, ...
    EEG.etc.alz.heartbeats.meanHr, replacedNote(replaced, opt.EventType));
end

% ======================================================================= %
function requireContinuous(input)
%REQUIRECONTINUOUS  The beats are found along one continuous recording.
    if ~ismatrix(input.data) || ~strcmpi(TransTools.FieldOr(input, 'DataFormat', 'CONTINUOUS'), 'CONTINUOUS')
        throw(MException('Alakazam:Heartbeats', ['I''m afraid Heartbeats works on the ' ...
            'continuous recording, and this dataset is %s. Would you run it before DefineBins?'], ...
            lower(char(string(TransTools.FieldOr(input, 'DataFormat', 'epoched'))))));
    end
end

function opt = normaliseOptions(options, input)
%NORMALISEOPTIONS  The defaults for anything a stored options struct lacks,
%   and a check that the channel is in this recording.
    opt.Channel       = char(string(TransTools.FieldOr(options, 'Channel', ecgChannelGuess(input))));
    opt.Polarity      = char(string(TransTools.FieldOr(options, 'Polarity', 'Automatic')));
    opt.MinIntervalMs = double(TransTools.FieldOr(options, 'MinIntervalMs', 300));
    opt.HighpassHz    = double(TransTools.FieldOr(options, 'HighpassHz', 0.5));
    opt.EventType     = strtrim(char(string(TransTools.FieldOr(options, 'EventType', 'heartbeat'))));
    labels = {input.chanlocs.labels};
    if isempty(opt.Channel) || ~any(strcmpi(labels, opt.Channel))
        throw(MException('Alakazam:Heartbeats', ['I''m afraid this recording has no channel ' ...
            '"%s" to find the heartbeats in. Would you choose its ECG channel?'], opt.Channel));
    end
    if ~any(strcmpi(opt.Polarity, {'Automatic', 'Upright', 'Inverted'}))
        throw(MException('Alakazam:Heartbeats', ['I''m afraid "%s" is not a polarity here: ' ...
            'Automatic, Upright or Inverted.'], opt.Polarity));
    end
    if ~(opt.MinIntervalMs > 0)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the shortest interval must be ' ...
            'more than 0 ms (spectHR''s is 300 ms, at most 200 beats a minute).']));
    end
    if ~(opt.HighpassHz >= 0) || opt.HighpassHz >= double(input.srate) / 2
        throw(MException('Alakazam:Heartbeats', ['I''m afraid a high-pass at %g Hz cannot ' ...
            'be used at %g Hz: it must lie below %g Hz (0 for none).'], opt.HighpassHz, ...
            double(input.srate), double(input.srate) / 2));
    end
    if isempty(opt.EventType)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the heartbeats need an event ' ...
            'code: "heartbeat", say.']));
    end
end

function label = ecgChannelGuess(input)
%ECGCHANNELGUESS  The first channel typed ECG, else the first whose label
%   says ECG or EKG (channelTypeFromLabel), else none.
    label = '';
    chanlocs = input.chanlocs;
    if isempty(chanlocs) || ~isfield(chanlocs, 'labels')
        return;
    end
    if isfield(chanlocs, 'type')
        typed = find(strcmpi(arrayfun(@(c) char(string(c.type)), chanlocs, 'UniformOutput', false), 'ECG'), 1);
        if ~isempty(typed)
            label = char(string(chanlocs(typed).labels));
            return;
        end
    end
    for k = 1:numel(chanlocs)
        if strcmpi(channelTypeFromLabel(char(string(chanlocs(k).labels))), 'ECG')
            label = char(string(chanlocs(k).labels));
            return;
        end
    end
end

function [polarity, decidedBy, skewness] = choosePolarity(ecg, srate, choice)
%CHOOSEPOLARITY  spectHR's check on all of the channel's valid samples, or
%   the user's word for it.
    skewness = NaN;
    switch lower(choice)
        case 'upright'
            polarity = 'normal';
            decidedBy = 'the user';
        case 'inverted'
            polarity = 'inverted';
            decidedBy = 'the user';
        otherwise
            valid = ecg(~isnan(ecg));
            [polarity, info] = Rpeak.polarity((0:numel(valid) - 1) / srate, valid);
            skewness = info.skewness;
            decidedBy = strrep(info.decidedBy, '_', ' ');
    end
end

function stretches = validStretches(EEG, ecg)
%VALIDSTRETCHES  [first last] of every stretch with no cut and no rejected
%   sample in it. A boundary at latency L joins sample floor(L) to the next.
    n = numel(ecg);
    breakAfter = false(1, n);
    if ~isempty(EEG.event) && isfield(EEG.event, 'type')
        types = arrayfun(@(e) char(string(e.type)), EEG.event, 'UniformOutput', false);
        cuts = [EEG.event(strcmpi(types, 'boundary')).latency];
        cuts = floor(double(cuts));
        cuts = cuts(cuts >= 1 & cuts < n);
        breakAfter(cuts) = true;
    end
    valid = ~isnan(ecg);
    starts = find(valid & [true, ~valid(1:end - 1) | breakAfter(1:end - 1)]);
    ends = find(valid & [~valid(2:end) | breakAfter(1:end - 1), true]);
    stretches = [starts(:), ends(:)];
end

function [EEG, replaced] = addBeatEvents(EEG, latency, ibi, eventType)
%ADDBEATEVENTS  One event per beat, in time order among the others, with
%   its IBI and HR. Heartbeats of the same code already there (from an
%   earlier run on this branch) are replaced.
    replaced = 0;
    events = EEG.event;
    if ~isempty(events) && isfield(events, 'type')
        types = arrayfun(@(e) char(string(e.type)), events, 'UniformOutput', false);
        old = strcmp(types, eventType);
        replaced = nnz(old);
        events = events(~old);
    end
    names = {'type', 'latency', 'duration', 'ibi', 'hr'};
    if ~isempty(events)
        names = union(fieldnames(events)', names, 'stable');
    end
    template = cell2struct(cell(numel(names), 1), names, 1);
    beats = repmat(template, 1, numel(latency));
    for k = 1:numel(latency)
        beats(k).type = eventType;
        beats(k).latency = latency(k);
        beats(k).duration = 0;
        beats(k).ibi = ibi(k);
        beats(k).hr = 60000 / ibi(k);
    end
    if isempty(events)
        events = beats;
    else
        for name = setdiff(names, fieldnames(events)')
            [events.(name{1})] = deal([]);
        end
        events = [reshape(orderfields(events, template), 1, []), beats];
    end
    [~, order] = sort([events.latency]);
    EEG.event = events(order);
end

function v = pluginVersion()
%PLUGINVERSION  This plugin's version, from its manifest.
    v = 'unknown';
    try
        manifest = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), 'Heartbeats.json')));
        v = manifest.Version;
    catch
    end
end

function word = upright2word(polarity)
    if strcmp(polarity, 'inverted')
        word = 'inverted';
    else
        word = 'upright';
    end
end

function text = replacedNote(replaced, eventType)
    text = '';
    if replaced > 0
        text = sprintf('; %d earlier "%s" event(s) replaced', replaced, eventType);
    end
end
