function [latency, ibiMs, polarityInfo, sync] = fromRecording(EEG, opt)
%FROMRECORDING  The R-peaks of an ECG recorded on another device, at their
%   places in the EEG, through the triggers both recordings received.
%
%   [LATENCY, IBIMS, POLARITYINFO, SYNC] = Sync.fromRecording(EEG, OPT):
%     - THE RECORDING is OPT.Recording (a path, or a name in the folder the
%       EEG was read from), or else the one file in that folder with the
%       EEG's own name and another format Alakazam reads (sub01.vhdr's
%       sub01.edf), as Eye tracking finds its .asc. Read as Alakazam reads
%       it.
%     - ITS ECG CHANNEL is OPT.RecordingChannel (a label), or else the one
%       typed or labelled ECG or EKG, or else its only channel.
%     - THE SHARED TRIGGERS are the events of the codes in OPT.SharedCodes
%       (blank: every code both recordings have), codes compared as
%       Sync.normaliseCode writes them. Sync.matchTriggers pairs them and
%       finds the line that carries the other recording's clock onto the
%       EEG's: the offset, and the drift between the two clocks.
%     - REFUSED when fewer than three triggers pair, when the triggers are
%       too regular to tell one offset from another, or when a pair lies
%       further from the line than OPT.MaxResidualMs (or one EEG sample, if
%       that is longer): the clocks would not be lined up well enough for
%       the beats.
%     - THE R-PEAKS are found in the other recording (Rpeak.beats) and
%       carried along the line into the EEG, each to its fraction of a
%       sample; a beat outside the EEG is counted and left out. IBIMS are
%       the other recording's intervals, on the EEG's clock.
%
%   Both recordings must be continuous (no cuts): one line carries one
%   clock onto the other. Run Heartbeats before SelectData.
%
%   SYNC records the file, the channel, the codes, the pairs and the
%   unpaired triggers on each side, the misfit (RMS and largest, ms), the
%   drift (ppm), the offset (s), and the beats found and left out.
%
%   See also SYNC.MATCHTRIGGERS, SYNC.FROMXDF, HEARTBEATS.
    if hasCuts(EEG)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid this recording has been cut ' ...
            '(it has boundary events), so one line cannot carry another recording''s clock ' ...
            'onto it. Run Heartbeats on the recording as it was read, before SelectData.']));
    end
    path = recordingFile(EEG, opt.Recording);
    other = readOther(path);
    if hasCuts(other)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the other recording, %s, has ' ...
            'cuts in it (boundary events), so one line cannot carry its clock onto this one.'], ...
            nameOf(path)));
    end
    [values, label] = channelOf(other, opt.RecordingChannel, path);

    [tHere, codesHere] = triggers(EEG, opt.EventType);
    [tThere, codesThere] = triggers(other, opt.EventType);
    wanted = sharedCodes(opt.SharedCodes);
    match = Sync.matchTriggers(tHere, codesHere, tThere, codesThere, wanted, 0.02);
    refuseUnlessGood(match, opt, double(EEG.srate), path);

    [beatTimes, ibiMs, polarityInfo, nStretches] = Rpeak.beats(values, double(other.srate), [], opt);
    latency = 1 + (match.a + match.b * beatTimes) * double(EEG.srate);
    inside = latency >= 1 & latency <= size(EEG.data, 2);
    latency = latency(inside);
    ibiMs = ibiMs(inside) * match.b;

    usedHere = sum(ismember(codesHere, match.shared));
    usedThere = sum(ismember(codesThere, match.shared));
    residualsMs = match.residuals * 1000;
    sync = struct( ...
        'route',          'recording', ...
        'file',           path, ...
        'channel',        label, ...
        'rate',           double(other.srate), ...
        'sharedCodes',    {match.shared}, ...
        'pairs',          numel(match.pairsHere), ...
        'unpairedHere',   usedHere - numel(match.pairsHere), ...
        'unpairedThere',  usedThere - numel(match.pairsThere), ...
        'residualRmsMs',  sqrt(mean(residualsMs .^ 2)), ...
        'residualMaxMs',  max(abs(residualsMs)), ...
        'driftPpm',       (match.b - 1) * 1e6, ...
        'offsetS',        match.a, ...
        'ambiguity',      match.ambiguity, ...
        'stretches',      nStretches, ...
        'beatsFound',     numel(beatTimes), ...
        'beatsUnplaced',  nnz(~inside));
end

% ======================================================================= %
function tf = hasCuts(EEG)
    tf = false;
    if isfield(EEG, 'event') && ~isempty(EEG.event) && isfield(EEG.event, 'type')
        types = arrayfun(@(e) char(string(e.type)), EEG.event, 'UniformOutput', false);
        tf = any(strcmpi(types, 'boundary'));
    end
end

function path = recordingFile(EEG, given)
%RECORDINGFILE  The other recording: as given (a path, or a name in the
%   EEG's own folder), or the one file beside the EEG with its name.
    raw = '';
    if isfield(EEG, 'etc') && isstruct(EEG.etc) && isfield(EEG.etc, 'alz') && isfield(EEG.etc.alz, 'rawFile')
        raw = char(EEG.etc.alz.rawFile);
    end
    [folder, base, ext] = fileparts(raw);
    given = strtrim(char(string(given)));
    if ~isempty(given)
        path = given;
        if ~isfile(path) && ~isempty(folder)
            path = fullfile(folder, given);
        end
        if ~isfile(path)
            throw(MException('Alakazam:Heartbeats', ['I''m afraid there is no recording %s ' ...
                '(nor in %s).'], given, folder));
        end
        return;
    end
    if isempty(raw)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid this recording does not say ' ...
            'which file it was read from, so the other recording cannot be found beside it. ' ...
            'Would you name it?']));
    end
    readable = setdiff(lower([rawFormats().extensions]), {'.erp'});   % an erpset is averaged
    found = dir(fullfile(folder, [base '.*']));
    names = {found(~[found.isdir]).name};
    [~, ~, exts] = cellfun(@fileparts, names, 'UniformOutput', false);
    keep = ismember(lower(exts), readable) & ~strcmpi(exts, ext);
    names = names(keep);
    if isempty(names)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid there is no other recording ' ...
            'named %s beside this one in %s. Would you name it?'], base, folder));
    end
    if numel(names) > 1
        throw(MException('Alakazam:Heartbeats', ['I''m afraid there are several recordings ' ...
            'named %s beside this one (%s). Would you name the one with the ECG?'], base, ...
            strjoin(names, ', ')));
    end
    path = fullfile(folder, names{1});
end

function other = readOther(path)
%READOTHER  The other recording, read as Alakazam reads its format.
    [folder, name, ext] = fileparts(path);
    switch lower(ext)
        case '.set'
            other = pop_loadset('filename', [name ext], 'filepath', folder);
        case '.vhdr'
            other = pop_loadbv(folder, [name ext]);
        case '.mat'
            other = load(path, 'EEG').EEG;
        otherwise
            other = readRecording(path);
    end
    if ~ismatrix(other.data) || size(other.data, 3) > 1
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the other recording, %s, is ' ...
            'epoched; it needs to be continuous.'], nameOf(path)));
    end
end

function [values, label] = channelOf(other, wanted, path)
%CHANNELOF  The other recording's ECG channel: by label, else typed or
%   labelled ECG or EKG, else its only channel.
    labels = arrayfun(@(c) char(string(c.labels)), other.chanlocs, 'UniformOutput', false);
    if isempty(labels)
        labels = arrayfun(@(k) sprintf('%d', k), 1:size(other.data, 1), 'UniformOutput', false);
    end
    wanted = strtrim(char(string(wanted)));
    k = [];
    if ~isempty(wanted)
        k = find(strcmpi(labels, wanted), 1);
    else
        if isfield(other.chanlocs, 'type')
            types = arrayfun(@(c) char(string(c.type)), other.chanlocs, 'UniformOutput', false);
            k = find(strcmpi(types, 'ECG'), 1);
        end
        if isempty(k)
            k = find(~cellfun(@isempty, regexpi(labels, '^(ECG|EKG)', 'once')), 1);
        end
        if isempty(k) && isscalar(labels)
            k = 1;
        end
    end
    if isempty(k)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the other recording, %s, has no ' ...
            'channel %s; its channels: %s. Would you name its ECG channel?'], nameOf(path), ...
            quotedOr(wanted, 'typed or labelled ECG'), strjoin(labels, ', ')));
    end
    values = double(other.data(k, :));
    label = labels{k};
end

function [times, codes] = triggers(EEG, heartbeatCode)
%TRIGGERS  A recording's events as triggers: times in seconds from its first
%   sample, codes as Sync.normaliseCode writes them; cuts and heartbeats left
%   out.
    times = zeros(1, 0);
    codes = cell(1, 0);
    if ~isfield(EEG, 'event') || isempty(EEG.event)
        return;
    end
    types = arrayfun(@(e) char(string(e.type)), EEG.event, 'UniformOutput', false);
    keep = ~strcmpi(types, 'boundary') & ~strcmp(types, heartbeatCode);
    times = ([EEG.event(keep).latency] - 1) / double(EEG.srate);
    codes = cellfun(@Sync.normaliseCode, types(keep), 'UniformOutput', false);
end

function codes = sharedCodes(text)
%SHAREDCODES  The codes asked for, as normalised codes; {} for all.
    parts = regexp(strtrim(char(string(text))), '[\s,;]+', 'split');
    parts = parts(~cellfun(@isempty, parts));
    codes = cellfun(@Sync.normaliseCode, parts, 'UniformOutput', false);
end

function refuseUnlessGood(match, opt, srate, path)
%REFUSEUNLESSGOOD  Say why the clocks cannot be lined up well enough.
    if ~isempty(match.problem)
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the two recordings could not be ' ...
            'lined up: %s (between this one and %s).'], match.problem, nameOf(path)));
    end
    pairs = numel(match.pairsHere);
    if pairs < 3
        throw(MException('Alakazam:Heartbeats', ['I''m afraid only %d shared trigger(s) could ' ...
            'be paired between this recording and %s; at least three are needed to line the ' ...
            'two clocks up and check it.'], pairs, nameOf(path)));
    end
    if match.ambiguity > 0.8
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the shared triggers (%s) are ' ...
            'too regular to tell how the two recordings line up: another offset fits them ' ...
            'nearly as well. Triggers with different codes, or at irregular times, settle it.'], ...
            strjoin(match.shared, ' ')));
    end
    limit = max(opt.MaxResidualMs, 1000 / srate);
    worst = max(abs(match.residuals)) * 1000;
    if worst > limit
        throw(MException('Alakazam:Heartbeats', ['I''m afraid the two recordings do not line up ' ...
            'well enough: after the best line through %d paired triggers, one is still %.2f ms ' ...
            'off, more than the %.2f ms allowed. A trigger one recording timed late, or clocks ' ...
            'that did not run steadily, would do that.'], pairs, worst, limit));
    end
end

function text = nameOf(path)
    [~, name, ext] = fileparts(path);
    text = [name ext];
end

function text = quotedOr(wanted, fallback)
    if isempty(wanted)
        text = fallback;
    else
        text = sprintf('"%s"', wanted);
    end
end
