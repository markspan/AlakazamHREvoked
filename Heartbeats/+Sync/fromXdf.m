function [latency, ibiMs, polarityInfo, sync] = fromXdf(EEG, opt)
%FROMXDF  The R-peaks of an ECG stream in the recording's own XDF file, at
%   their places in the EEG.
%
%   [LATENCY, IBIMS, POLARITYINFO, SYNC] = Sync.fromXdf(EEG, OPT) reads the
%   XDF file the recording was read from (EEG.etc.alz.rawFile) again, with
%   load_xdf and its defaults, as Alakazam's import does (pop_loadxdf):
%   every stream's time stamps on the recording computer's clock, each
%   stream's clock offsets fitted robustly with their drift, clock resets
%   handled, jitter removed, and lags the file declares corrected. So the
%   EEG stream's time stamps are the ones its samples were read with.
%
%   The ECG stream is OPT.Stream (its name), or else the first stream typed
%   ECG, or else one whose name says ECG or EKG; its channel OPT.StreamChannel
%   (a label), or else one labelled ECG or EKG, or else its first. Its
%   R-peaks are found on its own samples and time stamps (Rpeak.beats), and
%   each is placed in the EEG at the fraction of a sample its time falls at
%   between the EEG stream's time stamps. (pop_loadxdf puts a marker of
%   another stream on the nearest sample; an R-peak keeps its fraction.) A
%   beat outside the EEG stream, or in a gap in it, has no place and is
%   counted. IBIMS are the ECG stream's own intervals.
%
%   The recording must still have the EEG stream's time base: as many
%   samples, at its rate. Run Heartbeats before Resample or SelectData.
%
%   SYNC records the file, the streams, the channel, the rates, the beats
%   found and those left without a place, and how many clock offsets the
%   ECG stream's footer lists (LabRecorder lists there the measurements
%   load_xdf synchronises the stream on).
%
%   See also RPEAK.BEATS, SYNC.FROMRECORDING, HEARTBEATS.
    file = '';
    if isfield(EEG, 'etc') && isstruct(EEG.etc) && isfield(EEG.etc, 'alz') && isfield(EEG.etc.alz, 'rawFile')
        file = char(EEG.etc.alz.rawFile);
    end
    [~, ~, ext] = fileparts(file);
    if ~any(strcmpi(ext, {'.xdf', '.xdfz'}))
        throw(MException('Alakazam:Heartbeats', ['I''m afraid this recording was not read ' ...
            'from an XDF file (it was read from %s), so it has no streams to take the ECG ' ...
            'from. Choose its ECG channel, or another recording.'], displayName(file)));
    end
    if ~isfile(file)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the XDF file this recording ' ...
            'was read from, %s, is no longer there, so its ECG stream cannot be read.'], file));
    end
    if exist('load_xdf', 'file') ~= 2 && ~isempty(which('pop_loadxdf'))
        % xdfimport keeps load_xdf in its xdf folder, which its EEGLAB hook
        % (eegplugin_xdfimport) adds to the path; add it as the hook does.
        addpath(fullfile(fileparts(which('pop_loadxdf')), 'xdf'));
    end
    if exist('load_xdf', 'file') ~= 2
        throw(MException('Alakazam:Heartbeats', ['I''m afraid reading the XDF file needs ' ...
            'EEGLAB''s xdfimport plugin, which is not on the path. Alakazam installs it the ' ...
            'first time it reads an XDF file.']));
    end

    streams = load_xdf(file);
    [eegStream, eegIndex] = eegStreamOf(streams, EEG);
    eegTimes = double(eegStream.time_stamps);
    nominal = str2double(eegStream.info.nominal_srate);
    % Every sample, at about the stream's rate: pop_loadxdf gives the nominal
    % rate, File-IO (FieldTrip) the effective one, which differs by ppm.
    if numel(eegTimes) ~= size(EEG.data, 2) || abs(nominal - double(EEG.srate)) > 1e-3 * nominal
        throw(MException('Alakazam:Heartbeats', ['I''m afraid this recording no longer has ' ...
            'the time base of its XDF stream "%s" (%d samples at %g Hz there, %d at %g Hz ' ...
            'here), so the stream''s times cannot be placed in it. Run Heartbeats on the ' ...
            'recording as it was read, before Resample or SelectData.'], eegStream.info.name, ...
            numel(eegTimes), nominal, size(EEG.data, 2), double(EEG.srate)));
    end

    ecgStream = ecgStreamOf(streams, eegIndex, opt.Stream);
    [values, label] = channelOf(ecgStream, opt.StreamChannel);
    ecgTimes = double(ecgStream.time_stamps);
    [beatTimes, ibiMs, polarityInfo, nStretches] = Rpeak.beats(values, ecgTimes, [], opt);

    latency = interp1(eegTimes, 1:numel(eegTimes), beatTimes, 'linear', NaN);
    steps = diff(eegTimes);
    for g = find(steps > 2.5 * median(steps))
        latency(beatTimes > eegTimes(g) & beatTimes < eegTimes(g + 1)) = NaN;
    end
    placed = ~isnan(latency);
    latency = latency(placed);
    ibiMs = ibiMs(placed);

    sync = struct( ...
        'route',          'xdf', ...
        'file',           file, ...
        'eegStream',      char(string(eegStream.info.name)), ...
        'stream',         char(string(ecgStream.info.name)), ...
        'channel',        label, ...
        'nominalRate',    str2double(ecgStream.info.nominal_srate), ...
        'effectiveRate',  fieldOr(ecgStream.info, 'effective_srate', NaN), ...
        'clockOffsets',   clockOffsetCount(ecgStream), ...
        'stretches',      nStretches, ...
        'beatsFound',     numel(beatTimes), ...
        'beatsUnplaced',  nnz(~placed));
end

% ======================================================================= %
function [stream, index] = eegStreamOf(streams, EEG)
%EEGSTREAMOF  The stream the EEG was read from: by the name and type the
%   import recorded, else as pop_loadxdf chooses, the first typed EEG.
    name = '';
    if isfield(EEG.etc, 'info') && isstruct(EEG.etc.info) && isfield(EEG.etc.info, 'name')
        name = char(string(EEG.etc.info.name));
    end
    for index = 1:numel(streams)
        info = streams{index}.info;
        if strcmp(fieldOr(info, 'type', ''), 'EEG') && (isempty(name) || strcmp(fieldOr(info, 'name', ''), name))
            stream = streams{index};
            return;
        end
    end
    throw(MException('Alakazam:Heartbeats', ['I''m afraid the XDF file has no EEG stream%s ' ...
        'any more, which this recording was read from.'], quoted(name)));
end

function stream = ecgStreamOf(streams, eegIndex, wanted)
%ECGSTREAMOF  The ECG stream: by name, else typed ECG, else named so.
    wanted = strtrim(char(string(wanted)));
    candidates = {};
    for k = 1:numel(streams)
        s = streams{k};
        if k == eegIndex || ~isnumeric(s.time_series) || isempty(s.time_series) ...
                || ~(str2double(fieldOr(s.info, 'nominal_srate', '0')) > 0)
            continue;
        end
        candidates{end + 1} = s; %#ok<AGROW>
    end
    names = cellfun(@(s) char(string(fieldOr(s.info, 'name', ''))), candidates, 'UniformOutput', false);
    types = cellfun(@(s) char(string(fieldOr(s.info, 'type', ''))), candidates, 'UniformOutput', false);
    if ~isempty(wanted)
        hit = find(strcmp(names, wanted), 1);
    else
        hit = find(strcmpi(types, 'ECG'), 1);
        if isempty(hit)
            hit = find(~cellfun(@isempty, regexpi(names, 'ECG|EKG', 'once')), 1);
        end
    end
    if isempty(hit)
        listed = strjoin(cellfun(@(n, t) sprintf('"%s" (%s)', n, t), names, types, ...
            'UniformOutput', false), ', ');
        if isempty(listed)
            listed = 'none';
        end
        if isempty(wanted)
            what = 'no stream typed ECG or named so';
        else
            what = sprintf('no stream "%s"', wanted);
        end
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the XDF file has %s. Its ' ...
            'signal streams besides the EEG: %s. Would you name the ECG stream?'], what, listed));
    end
    stream = candidates{hit};
end

function [values, label] = channelOf(stream, wanted)
%CHANNELOF  The ECG channel of a stream: by label (or number), else the one
%   labelled ECG or EKG, else the first.
    labels = streamLabels(stream);
    wanted = strtrim(char(string(wanted)));
    if ~isempty(wanted)
        k = find(strcmpi(labels, wanted), 1);
        if isempty(k) && ~isnan(str2double(wanted)) && str2double(wanted) <= numel(labels)
            k = str2double(wanted);
        end
        if isempty(k)
            throw(MException('Alakazam:Heartbeats', ['I''m afraid the stream "%s" has no ' ...
                'channel "%s"; its channels: %s.'], char(string(stream.info.name)), wanted, ...
                strjoin(labels, ', ')));
        end
    else
        k = find(~cellfun(@isempty, regexpi(labels, '^(ECG|EKG)', 'once')), 1);
        if isempty(k)
            k = 1;
        end
    end
    values = double(stream.time_series(k, :));
    label = labels{k};
end

function labels = streamLabels(stream)
%STREAMLABELS  A stream's channel labels, or numbers where it has none.
    n = size(stream.time_series, 1);
    labels = arrayfun(@(k) sprintf('%d', k), 1:n, 'UniformOutput', false);
    try
        channels = stream.info.desc.channels.channel;
        for k = 1:min(n, numel(channels))
            if iscell(channels)
                c = channels{k};
            else
                c = channels(k);
            end
            if isfield(c, 'label') && ~isempty(c.label)
                labels{k} = char(string(c.label));
            end
        end
    catch
    end
end

function n = clockOffsetCount(stream)
    n = 0;
    try
        n = numel(stream.info.clock_offsets.offset);
    catch
    end
end

function v = fieldOr(s, name, default)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end

function text = displayName(file)
    if isempty(file)
        text = 'a file it does not name';
    else
        [~, name, ext] = fileparts(file);
        text = [name ext];
    end
end

function text = quoted(name)
    text = '';
    if ~isempty(name)
        text = sprintf(' "%s"', name);
    end
end
