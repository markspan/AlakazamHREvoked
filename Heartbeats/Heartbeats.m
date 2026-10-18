function [EEG, options] = Heartbeats(input, varargin)
%% HEARTBEATS  Every heartbeat as an event, from the R-peaks of an ECG: the
%   first step to heartbeat-evoked potentials.
%
%   THE ECG can come from three places (Source):
%     'channel'    a channel of this recording;
%     'xdf'        a stream of the XDF file this recording was read from,
%                  its times on the recording computer's clock as load_xdf
%                  puts them (Sync.fromXdf);
%     'recording'  a recording made by another device that received the
%                  same triggers, its clock lined up with this one's on them,
%                  offset and drift (Sync.fromRecording).
%   Only the beats come over, never the ECG itself: each R-peak is found in
%   the ECG's own samples and placed in the EEG at the fraction of a sample
%   its time falls at.
%
%   The R-peaks are found exactly as spectHR finds them, ported to MATLAB
%   (the +Rpeak package, held to spectHR's own results by
%   SpectHRParityTest):
%     1. the ECG's polarity, from the skewness of its 5 to 20 Hz band (an
%        inverted ECG is turned upright for detection, not in the data);
%     2. a zero-phase 0.5 Hz high-pass, order 4, against baseline wander;
%     3. the local maxima at least the median plus 1.5 standard deviations,
%        no two closer than the shortest interval (300 ms), each timed
%        between samples, at the apex of a symmetric triangle through it and
%        its two neighbours (Rpeak.detect).
%   The data are not changed. Each R-peak becomes an event of the code
%   chosen ('heartbeat' by default), at its sub-sample latency, with two
%   fields: IBI, the interval to the next beat in ms (as spectHR's IBI
%   series), and HR, the heart rate it makes in beats per minute; NaN for
%   the last beat before a cut, a gap or the end.
%
%   THE CHECK (Check, on by default): the EEG averaged at the beats shows
%   the heart's own electrical field, which peaks with the R-wave; where it
%   peaks relative to the beats is the timing's error from end to end
%   (Sync.heartField). It is recorded and, with CorrectLag, the beats are
%   moved by it.
%
%   THE BINS come next, from DefineBins, which can choose heartbeats by what
%   happens around them, e.g. the heartbeats in the second after a face:
%       bin 1 "HEP after faces" : "heartbeat" any(1) within [-1000, 0] ms
%   and Average then gives the heartbeat-evoked potential of each bin.
%
%   CUTS AND REJECTED STRETCHES. An ECG channel cut by 'boundary' events, or
%   with rejected (NaN) samples, is a set of separate stretches: each is
%   filtered and searched on its own, as spectHR would a recording of its
%   own, so no beat is found across a join. The polarity is decided once,
%   on all of the ECG's valid samples.
%
%   THE ECG's ELECTRICAL FIELD in the EEG (the cardiac field artefact) is
%   not removed here: remove the heart components with ICA after.
%
%   WHAT IT RECORDS, in EEG.etc.alz.heartbeats: the plugin version, the
%   source and its synchronisation (sync), the channel, the polarity and
%   what decided it, the settings, the beats, the mean heart rate and the
%   check.
%
%   Signature (Alakazam transformation contract):
%     [EEG, options] = Heartbeats(input)        % interactive dialog
%     [EEG, options] = Heartbeats(input, opts)  % replay a stored options struct
%
%   See also RPEAK.BEATS, SYNC.FROMXDF, SYNC.FROMRECORDING, SYNC.HEARTFIELD,
%   DEFINEBINS.
[opts, interactive] = TransTools.InitGuard(nargin, 'Alakazam:Heartbeats', varargin{:});

requireContinuous(input);

if interactive
    stored = TransformSettings.get('Heartbeats');
    if isempty(stored) || ~isstruct(stored)
        stored = struct();
    end
    d = @(f, v) TransTools.FieldOr(stored, f, v);
    from = @(source) @(v) strcmp(v.Source, source);
    options = TransformOptionsDialog( ...
        'title', 'Heartbeats', ...
        'Description', ['Every R-peak of the ECG becomes an event, found as spectHR finds it. ' ...
            'The data are not changed. Choose the heartbeats per condition with DefineBins ' ...
            'afterwards, and Average them.'], ...
        {'ECG from'; 'Source'}, DialogFields.Choice({'A channel of this recording', ...
            'A stream in this recording''s XDF file', 'Another recording, on shared triggers'}, ...
            d('Source', 'channel'), 'Values', {'channel', 'xdf', 'recording'}), ...
        {'ECG channel'; 'Channel'}, DialogFields.Channels(input, d('Channel', ecgChannelGuess(input)), ...
            'Multiple', false, 'Required', true, 'QuickSelect', cell(0, 2), ...
            'EnabledWhen', from('channel')), ...
        {'Stream (blank: typed ECG)'; 'Stream'}, DialogFields.Text(d('Stream', ''), ...
            'EnabledWhen', from('xdf')), ...
        {'Its channel (blank: ECG)'; 'StreamChannel'}, DialogFields.Text(d('StreamChannel', ''), ...
            'EnabledWhen', from('xdf')), ...
        {'Recording (blank: same name)'; 'Recording'}, DialogFields.Text(d('Recording', ''), ...
            'EnabledWhen', from('recording'), 'Tooltip', ['The other device''s recording; blank: ' ...
            'the file beside this one with its name and another format']), ...
        {'Its ECG channel (blank: ECG)'; 'RecordingChannel'}, DialogFields.Text(d('RecordingChannel', ''), ...
            'EnabledWhen', from('recording')), ...
        {'Shared triggers (blank: all)'; 'SharedCodes'}, DialogFields.Text(d('SharedCodes', ''), ...
            'EnabledWhen', from('recording'), 'Tooltip', ['The codes of the triggers both ' ...
            'recordings received, e.g. 1 2 3; blank: every code both have']), ...
        {'Largest trigger misfit (ms)'; 'MaxResidualMs'}, DialogFields.Number(d('MaxResidualMs', 2), ...
            'Limits', [0.01 Inf], 'EnabledWhen', from('recording')), ...
        'separator', 'Detection, as spectHR detects:', ...
        {'ECG polarity'; 'Polarity'}, DialogFields.Choice({'Automatic', 'Upright', 'Inverted'}, ...
            d('Polarity', 'Automatic')), ...
        {'Shortest interval (ms)'; 'MinIntervalMs'}, DialogFields.Number(d('MinIntervalMs', 300), ...
            'Limits', [1 Inf]), ...
        {'High-pass before detection (Hz)'; 'HighpassHz'}, DialogFields.Number(d('HighpassHz', 0.5), ...
            'Limits', [0 Inf]), ...
        {'Event code'; 'EventType'}, DialogFields.Text(d('EventType', 'heartbeat')), ...
        'separator', 'Timing:', ...
        {'Check on the heart''s field in the EEG'; 'Check'}, DialogFields.Checkbox(d('Check', true)), ...
        {'Move the beats by the lag found'; 'CorrectLag'}, DialogFields.Checkbox(d('CorrectLag', false), ...
            'EnabledWhen', @(v) v.Check));
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
srate = double(EEG.srate);
ecgChannel = [];
switch opt.Source
    case 'channel'
        ecgChannel = find(strcmpi({EEG.chanlocs.labels}, opt.Channel), 1);
        [beatTimes, ibi, polarity, nStretches] = Rpeak.beats(double(EEG.data(ecgChannel, :)), ...
            srate, cutsOf(EEG), opt);
        latency = 1 + beatTimes * srate;
        sync = struct('route', 'channel', 'channel', opt.Channel, 'stretches', nStretches);
    case 'xdf'
        [latency, ibi, polarity, sync] = Sync.fromXdf(EEG, opt);
    case 'recording'
        [latency, ibi, polarity, sync] = Sync.fromRecording(EEG, opt);
end
if isempty(latency)
    throw(MException('Alakazam:Heartbeats', ['I''m afraid no heartbeat was found in %s. ' ...
        'Is it the ECG, and is the ECG in it?'], sync.channel));
end

check = struct('seen', false, 'lagMs', NaN, 'ratio', NaN, 'beats', 0, 'channels', 0, ...
    'corrected', false, 'movedOut', 0);
if opt.Check
    check = Sync.heartField(EEG, latency, ecgChannel);
    check.corrected = opt.CorrectLag && check.seen;
    check.movedOut = 0;
    if check.corrected
        latency = latency + check.lagMs / 1000 * srate;
        inside = latency >= 1 & latency <= size(EEG.data, 2);   % a beat moved past either end
        check.movedOut = nnz(~inside);
        latency = latency(inside);
        ibi = ibi(inside);
    end
end

[EEG, replaced] = addBeatEvents(EEG, latency, ibi, opt.EventType);
hr = 60000 ./ ibi;
EEG.etc.alz.heartbeats = struct( ...
    'version',        pluginVersion(), ...
    'source',         opt.Source, ...
    'channel',        sync.channel, ...
    'polarity',       polarity.polarity, ...
    'polarityBy',     polarity.decidedBy, ...
    'skewness',       polarity.skewness, ...
    'minIntervalMs',  opt.MinIntervalMs, ...
    'highpassHz',     opt.HighpassHz, ...
    'eventType',      opt.EventType, ...
    'beats',          numel(latency), ...
    'stretches',      sync.stretches, ...
    'meanHr',         mean(hr(isfinite(hr))), ...
    'replaced',       replaced, ...
    'sync',           sync, ...
    'check',          check, ...
    'options',        options);

fprintf('Heartbeats: %d heartbeats in %s (%s ECG, by %s), %.1f beats a minute on average%s.\n', ...
    numel(latency), sourceName(sync), upright2word(polarity.polarity), polarity.decidedBy, ...
    EEG.etc.alz.heartbeats.meanHr, replacedNote(replaced, opt.EventType));
reportSync(sync);
reportCheck(check);
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
%NORMALISEOPTIONS  The defaults for anything a stored options struct lacks
%   (one stored before a setting existed replays as it ran: from a channel,
%   with the check), and a check that what is asked for can be done.
    opt.Source           = lower(char(string(TransTools.FieldOr(options, 'Source', 'channel'))));
    opt.Channel          = char(string(TransTools.FieldOr(options, 'Channel', ecgChannelGuess(input))));
    opt.Stream           = char(string(TransTools.FieldOr(options, 'Stream', '')));
    opt.StreamChannel    = char(string(TransTools.FieldOr(options, 'StreamChannel', '')));
    opt.Recording        = char(string(TransTools.FieldOr(options, 'Recording', '')));
    opt.RecordingChannel = char(string(TransTools.FieldOr(options, 'RecordingChannel', '')));
    opt.SharedCodes      = char(string(TransTools.FieldOr(options, 'SharedCodes', '')));
    opt.MaxResidualMs    = double(TransTools.FieldOr(options, 'MaxResidualMs', 2));
    opt.Polarity         = char(string(TransTools.FieldOr(options, 'Polarity', 'Automatic')));
    opt.MinIntervalMs    = double(TransTools.FieldOr(options, 'MinIntervalMs', 300));
    opt.HighpassHz       = double(TransTools.FieldOr(options, 'HighpassHz', 0.5));
    opt.EventType        = strtrim(char(string(TransTools.FieldOr(options, 'EventType', 'heartbeat'))));
    opt.Check            = logical(TransTools.FieldOr(options, 'Check', true));
    opt.CorrectLag       = logical(TransTools.FieldOr(options, 'CorrectLag', false));
    if ~any(strcmp(opt.Source, {'channel', 'xdf', 'recording'}))
        throw(MException('Alakazam:Heartbeats', ['I''m afraid "%s" is not where an ECG can ' ...
            'come from here: channel, xdf or recording.'], opt.Source));
    end
    labels = {input.chanlocs.labels};
    if strcmp(opt.Source, 'channel') && (isempty(opt.Channel) || ~any(strcmpi(labels, opt.Channel)))
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
    if ~(opt.HighpassHz >= 0) || (strcmp(opt.Source, 'channel') && opt.HighpassHz >= double(input.srate) / 2)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid a high-pass at %g Hz cannot ' ...
            'be used at %g Hz: it must lie below %g Hz (0 for none).'], opt.HighpassHz, ...
            double(input.srate), double(input.srate) / 2));
    end
    if ~(opt.MaxResidualMs > 0)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the largest trigger misfit must ' ...
            'be more than 0 ms.']));
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

function breakAfter = cutsOf(EEG)
%CUTSOF  Per sample, whether a cut follows it: a boundary at latency L joins
%   sample floor(L) to the next.
    n = size(EEG.data, 2);
    breakAfter = false(1, n);
    if ~isempty(EEG.event) && isfield(EEG.event, 'type')
        types = arrayfun(@(e) char(string(e.type)), EEG.event, 'UniformOutput', false);
        cuts = floor(double([EEG.event(strcmpi(types, 'boundary')).latency]));
        breakAfter(cuts(cuts >= 1 & cuts < n)) = true;
    end
end

function text = sourceName(sync)
%SOURCENAME  Where the ECG came from, for the console.
    switch sync.route
        case 'xdf'
            text = sprintf('the XDF stream "%s" (%s)', sync.stream, sync.channel);
        case 'recording'
            [~, name, ext] = fileparts(sync.file);
            text = sprintf('%s (%s)', [name ext], sync.channel);
        otherwise
            text = sync.channel;
    end
end

function reportSync(sync)
%REPORTSYNC  How the other clock was lined up, for the console.
    switch sync.route
        case 'xdf'
            fprintf(['Heartbeats: the stream''s footer lists %d clock offset(s); %d of %d beats ' ...
                'outside the EEG stream.\n'], sync.clockOffsets, sync.beatsUnplaced, sync.beatsFound);
        case 'recording'
            fprintf(['Heartbeats: %d shared triggers paired (%d unpaired here, %d there); misfit ' ...
                '%.2f ms RMS, %.2f ms at most; drift %.1f ppm; %d of %d beats outside the EEG.\n'], ...
                sync.pairs, sync.unpairedHere, sync.unpairedThere, sync.residualRmsMs, ...
                sync.residualMaxMs, sync.driftPpm, sync.beatsUnplaced, sync.beatsFound);
    end
end

function reportCheck(check)
%REPORTCHECK  The check on the heart's field, for the console.
    if check.beats == 0 && check.channels == 0
        return;
    end
    if ~check.seen
        fprintf(['Heartbeats: the heart''s field is not clear enough in the EEG (%d beats, %d ' ...
            'channels) to check the timing on.\n'], check.beats, check.channels);
        return;
    end
    fprintf('Heartbeats: the heart''s field in the EEG peaks %+.1f ms from the beats', check.lagMs);
    if check.corrected
        fprintf('; the beats were moved onto it (%d moved past an end, left out).\n', check.movedOut);
    elseif abs(check.lagMs) > 10
        fprintf('. WARNING: more than 10 ms; the ECG may lag the EEG (CorrectLag moves the beats).\n');
    else
        fprintf('.\n');
    end
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
