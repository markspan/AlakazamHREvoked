classdef SyncTest < matlab.unittest.TestCase
%SYNCTEST  Heartbeats from an ECG recorded apart from the EEG: the beats
%   carried onto the EEG's clock, and the check on the heart's own field.
%
%   Every recording here is made up, its heartbeats and clocks known:
%     - an EEG at 500 Hz whose scalp channels carry the heart's field (the
%       ECG, weighted per channel, as volume conduction puts it there) and
%       triggers at irregular times;
%     - an XDF file, written by the test itself, with that EEG and an ECG as
%       two streams, each on its own clock and with jittered time stamps, as
%       Lab Streaming Layer records them (route 'xdf');
%     - a second recording, saved as an EEGLAB dataset beside the EEG, with
%       the ECG at 1000 Hz on a clock offset by seconds and drifting by
%       40 ppm, some triggers missing and some extra (route 'recording').
%   The beats must land where they are, the check must find the field at
%   the beats, find a device's delay when there is one, and move the beats
%   by it when asked; and what cannot be lined up must be refused.
%
%   Run with: runtests('tests/SyncTest.m').
%
%   See also SYNC.MATCHTRIGGERS, SYNC.FROMXDF, SYNC.FROMRECORDING,
%   SYNC.HEARTFIELD, HEARTBEATSTEST.

    properties (Constant)
        Srate = 500
        Seconds = 90
        Detection = struct('Polarity', 'Automatic', 'MinIntervalMs', 300, 'HighpassHz', 0.5)
    end

    properties
        Folder
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            plugin = fileparts(fileparts(mfilename('fullpath')));
            root = getenv('ALAKAZAM_ROOT');   % Alakazam's source, which the plugin runs in
            if isempty(root)
                root = fullfile(fileparts(plugin), 'Alakazam');
            end
            testCase.assumeTrue(isfolder(fullfile(root, 'src')), ...
                'Alakazam''s source is not found (set ALAKAZAM_ROOT).');
            for p = {fullfile(root, 'src'), fullfile(root, 'src', 'Transformations'), ...
                     fullfile(root, 'src', 'Support'), fullfile(root, 'src', 'IO'), ...
                     fullfile(root, 'src', 'Dialogs'), fullfile(plugin, 'Heartbeats')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(p{1}));
            end
        end
    end

    methods (TestMethodSetup)
        function ownFolder(testCase)
            testCase.Folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        % ---- the pieces -------------------------------------------------
        function triggersPairDespiteMissingAndExtraOnes(testCase)
        %TRIGGERSPAIRDESPITEMISSINGANDEXTRAONES  An offset of 12.3 s and a
        %   drift of 40 ppm found from triggers of three codes over an hour
        %   (the drift then 144 ms, far more than the pairing's tolerance),
        %   three missing there and two extra, one of them a spurious trigger
        %   0.7 s from a lost one, of its code; each timed to within 0.1 ms.
            stream = RandStream('mt19937ar', 'Seed', 3);
            tHere = cumsum(10 + 40 * rand(stream, 1, 120));
            codesHere = arrayfun(@(k) sprintf('%d', mod(k, 3) + 1), 1:120, 'UniformOutput', false);
            a = 12.3;
            b = 1 + 40e-6;
            tThere = (tHere - a) / b + 1e-4 * (rand(stream, 1, 120) - 0.5);
            codesThere = codesHere;
            spurious = tThere(60) + 0.7;
            missing = [5 60 99];
            tThere(missing) = [];
            codesThere(missing) = [];
            tThere = [tThere, spurious, 150];
            codesThere = [codesThere, codesHere(60), {'2'}];

            m = Sync.matchTriggers(tHere, codesHere, tThere, codesThere, {}, 0.02);
            testCase.assertEmpty(m.problem);
            testCase.verifyEqual(numel(m.pairsHere), 117);
            testCase.verifyEqual(m.a, a, 'AbsTol', 1e-3);
            testCase.verifyEqual((m.b - 1) * 1e6, 40, 'AbsTol', 1);
            testCase.verifyLessThan(max(abs(m.residuals)), 1e-4);
            testCase.verifyLessThan(m.ambiguity, 0.5);
        end

        function fewTriggersFarApartLineUpToo(testCase)
        %FEWTRIGGERSFARAPARTLINEUPTOO  Only a block's onset shared, twelve
        %   over two hours, the clocks drifting 150 ppm (over a second in
        %   all): too far apart for the pairing to follow the drift from one
        %   to the next, so the vote over rate ratios has to find it.
            stream = RandStream('mt19937ar', 'Seed', 4);
            tHere = cumsum(500 + 200 * rand(stream, 1, 12));
            b = 1 + 150e-6;
            tThere = (tHere - 3.1) / b + 1e-4 * (rand(stream, 1, 12) - 0.5);
            m = Sync.matchTriggers(tHere, repmat({'block'}, 1, 12), tThere, repmat({'block'}, 1, 12), {}, 0.02);
            testCase.assertEmpty(m.problem);
            testCase.verifyEqual(numel(m.pairsHere), 12);
            testCase.verifyEqual((m.b - 1) * 1e6, 150, 'AbsTol', 0.1);
        end

        function tooRegularTriggersAreAmbiguous(testCase)
        %TOOREGULARTRIGGERSAREAMBIGUOUS  One code every 2 s exactly: an offset
        %   by whole periods fits nearly as well, which the vote says.
            tHere = 2 * (1:100);
            m = Sync.matchTriggers(tHere, repmat({'1'}, 1, 100), tHere - 7.25, repmat({'1'}, 1, 100), {}, 0.02);
            testCase.verifyGreaterThan(m.ambiguity, 0.9);
        end

        function anEcgStreamWithAGapIsSearchedInTwoStretches(testCase)
        %ANECGSTREAMWITHAGAPISSEARCHEDINTWOSTRETCHES  Four beats of a 1000 Hz
        %   ECG stream lost (its time stamps jump): two stretches, every beat
        %   either side found, no interval across the gap.
            stream = RandStream('mt19937ar', 'Seed', 8);
            t = 100 + (0:59999) / 1000;
            beatTimes = SyncTest.beats(100.6, 159, stream);
            ecg = SyncTest.ecgAt(t, beatTimes) + 10 * randn(stream, size(t));
            lost = t > beatTimes(30) + 0.4 & t < beatTimes(34) + 0.4;
            [found, ibi, ~, stretches] = Rpeak.beats(ecg(~lost), t(~lost), [], SyncTest.Detection);
            testCase.verifyEqual(stretches, 2);
            testCase.verifyEqual(found, beatTimes([1:30, 35:end]), 'AbsTol', 0.002);
            testCase.verifyEqual(find(isnan(ibi)), [30, numel(found)], 'No interval across the gap.');
        end

        function aRobustLineIgnoresAFewWildPoints(testCase)
            x = 1:50;
            y = 3 + 2 * x;
            y([7 31]) = y([7 31]) + [40, -25];
            [a, b, r] = Sync.robustLine(x, y);
            testCase.verifyEqual([a b], [3 2], 'AbsTol', 1e-6);
            testCase.verifyEqual(abs(r([7 31])), [40 25], 'AbsTol', 1e-5);
        end

        function codesAreReadAsDefineBinsReadsThem(testCase)
            testCase.verifyEqual(Sync.normaliseCode('S  1'), '1');
            testCase.verifyEqual(Sync.normaliseCode('007'), '7');
            testCase.verifyEqual(Sync.normaliseCode('R 12'), 'R12');
            testCase.verifyEqual(Sync.normaliseCode(5), '5');
            testCase.verifyEqual(Sync.normaliseCode('sync'), 'sync');
        end

        function theCheckFindsTheFieldAtTheBeatsAndALagWhereThereIsOne(testCase)
            [EEG, truth] = SyncTest.eeg();
            check = Sync.heartField(EEG, truth, []);
            testCase.verifyTrue(check.seen);
            testCase.verifyEqual(check.lagMs, 0, 'AbsTol', 1.5);
            late = Sync.heartField(EEG, truth + 0.025 * SyncTest.Srate, []);
            testCase.verifyEqual(late.lagMs, -25, 'AbsTol', 0.5, ...
                'Beats 25 ms late: the field 25 ms before them, between samples.');

            EEG.data(1:4, :) = randn(RandStream('mt19937ar', 'Seed', 1), 4, size(EEG.data, 2));
            testCase.verifyFalse(Sync.heartField(EEG, truth, []).seen, 'No field, no conclusion.');
        end

        function anEarlyEcgChannelIsMovedOntoTheField(testCase)
        %ANEARLYECGCHANNELISMOVEDONTOTHEFIELD  An ECG channel 40 ms ahead of
        %   the EEG channels (a beat just after the recording's end shows in
        %   it 20 ms before the end): the check finds the lag, CorrectLag
        %   moves the beats onto the field, and the one moved past the last
        %   sample is left out.
            [EEG, truth] = SyncTest.eeg();
            t = (0:size(EEG.data, 2) - 1) / SyncTest.Srate;
            beatTimes = [(truth - 1) / SyncTest.Srate, SyncTest.Seconds + 0.02];
            ecg = SyncTest.ecgAt(t, beatTimes - 0.040) + 10 * randn(RandStream('mt19937ar', 'Seed', 2), size(t));
            EEG.data(5, :) = single(ecg);
            EEG.nbchan = 5;
            EEG.chanlocs(5) = struct('labels', 'ECG', 'type', 'ECG');

            out = Heartbeats(EEG, struct('Channel', 'ECG', 'CorrectLag', true));
            check = out.etc.alz.heartbeats.check;
            testCase.verifyEqual(check.lagMs, 40, 'AbsTol', 0.5);
            testCase.verifyTrue(check.corrected);
            testCase.verifyEqual(check.movedOut, 1, 'The beat 20 ms before the end, moved 20 ms past it.');
            testCase.verifyEqual(SyncTest.lateBy(out, truth), zeros(1, numel(truth)), 'AbsTol', 2, ...
                'Each beat within the detector''s precision of the true one.');
        end

        % ---- route 'xdf' ------------------------------------------------
        function anEcgStreamsBeatsLandWhereTheyAreInTheEeg(testCase)
        %ANECGSTREAMSBEATSLANDWHERETHEYAREINTHEEEG  Two streams on clocks
        %   1.25 s and -3.7 s off the recorder's, jittered time stamps; the
        %   ECG at 1000 Hz starting before the EEG and ending after it. Each
        %   beat the detector finds lands where the true clocks put it, to
        %   within a hundredth of a sample; and so within the detector's own
        %   precision (on this noisy ECG, about a ms) of the true beat.
            testCase.requireXdf();
            [EEG, truthLatency, detected] = testCase.xdfRecording(1000, 0);
            xdfFolder = fileparts(which('load_xdf'));
            rmpath(xdfFolder);   % as where xdfimport's EEGLAB hook has not run
            testCase.addTeardown(@addpath, xdfFolder);
            out = Heartbeats(EEG, struct('Source', 'xdf'));

            beats = out.event(strcmp({out.event.type}, 'heartbeat'));
            inside = detected >= 1 & detected <= size(EEG.data, 2);
            testCase.assertNumElements(beats, nnz(inside));
            testCase.verifyLessThan(max(abs([beats.latency] - detected(inside))), 0.01, ...
                'Carried onto the EEG as the true clocks carry them.');
            testCase.verifyLessThan(max(abs([beats.latency] - truthLatency(inside))), 1, ...
                'Within a sample (2 ms) of the true beats.');
            ibi = [beats.ibi];
            testCase.verifyEqual(ibi(1:end - 1), diff(detected(inside)) / SyncTest.Srate * 1000, ...
                'AbsTol', 0.02, 'The intervals, on the ECG''s own samples.');
            sync = out.etc.alz.heartbeats.sync;
            testCase.verifyEqual(sync.route, 'xdf');
            testCase.verifyEqual(sync.stream, 'Chest ECG');
            testCase.verifyEqual(sync.clockOffsets, 29);
            testCase.verifyEqual(sync.beatsUnplaced, nnz(~inside), 'The beats outside the EEG, counted.');
            check = out.etc.alz.heartbeats.check;
            testCase.verifyTrue(check.seen);
            testCase.verifyEqual(check.lagMs, 0, 'AbsTol', 1.5);
            testCase.verifyEqual(out.data, EEG.data, 'The data are not changed.');
        end

        function aLateDeviceIsFoundAndCorrectedOnRequest(testCase)
        %ALATEDEVICEISFOUNDANDCORRECTEDONREQUEST  A 130 Hz chest strap whose
        %   samples are stamped 40 ms late: the check finds the field 40 ms
        %   before the beats, and with CorrectLag the beats move onto it.
            testCase.requireXdf();
            [EEG, truthLatency] = testCase.xdfRecording(130, 0.040);

            asIs = Heartbeats(EEG, struct('Source', 'xdf'));
            testCase.verifyEqual(asIs.etc.alz.heartbeats.check.lagMs, -40, 'AbsTol', 3);
            testCase.verifyEqual(median(SyncTest.lateBy(asIs, truthLatency)), 40, 'AbsTol', 2, ...
                'Placed as stamped: 40 ms late.');

            moved = Heartbeats(EEG, struct('Source', 'xdf', 'CorrectLag', true));
            testCase.verifyTrue(moved.etc.alz.heartbeats.check.corrected);
            testCase.verifyEqual(median(SyncTest.lateBy(moved, truthLatency)), 0, 'AbsTol', 2, ...
                'Moved onto the heart''s field.');
            latency = [moved.event(strcmp({moved.event.type}, 'heartbeat')).latency];
            testCase.verifyTrue(all(latency >= 1 & latency <= size(EEG.data, 2)), ...
                'None moved out of the recording.');
        end

        function aBeatInAGapOfTheEegStreamHasNoPlace(testCase)
        %ABEATINAGAPOFTHEEEGSTREAMHASNOPLACE  Two seconds of the EEG stream
        %   lost (a dropped connection; the import joins the samples either
        %   side): the beats in them are left out and counted, the others
        %   placed by the time stamps either side.
            testCase.requireXdf();
            [EEG, ~, detected, inGap] = testCase.xdfRecording(1000, 0, 20001:21000);
            out = Heartbeats(EEG, struct('Source', 'xdf'));
            beats = out.event(strcmp({out.event.type}, 'heartbeat'));
            placed = isfinite(detected);
            testCase.assertNumElements(beats, nnz(placed));
            testCase.verifyLessThan(max(abs([beats.latency] - detected(placed))), 0.01);
            testCase.assertGreaterThan(nnz(inGap), 0, 'Some beats fall in the gap.');
            testCase.verifyEqual(out.etc.alz.heartbeats.sync.beatsUnplaced, nnz(~placed));
        end

        function anXdfReadByFileIoWorksToo(testCase)
        %ANXDFREADBYFILEIOWORKSTOO  Read through File-IO, the EEG stream has
        %   its effective rate (here 30 ppm off the nominal) and no stream
        %   information; its samples and their time stamps are the same.
            testCase.requireXdf();
            [EEG, ~, detected] = testCase.xdfRecording(1000, 0);
            EEG.srate = SyncTest.Srate * (1 + 30e-6);
            EEG.etc = rmfield(EEG.etc, 'info');
            out = Heartbeats(EEG, struct('Source', 'xdf'));
            inside = detected >= 1 & detected <= size(EEG.data, 2);
            testCase.verifyNumElements(out.event(strcmp({out.event.type}, 'heartbeat')), nnz(inside));
        end

        function whatAnXdfCannotGiveIsRefused(testCase)
            testCase.requireXdf();
            EEG = testCase.xdfRecording(1000, 0);
            testCase.verifySubstring(SyncTest.refusal(@() Heartbeats(EEG, struct('Source', 'xdf', ...
                'Stream', 'Polar'))), 'no stream "Polar"');
            cut = EEG;
            cut.data = cut.data(:, 1:end - 10);
            cut.pnts = size(cut.data, 2);
            testCase.verifySubstring(SyncTest.refusal(@() Heartbeats(cut, struct('Source', 'xdf'))), ...
                'no longer has the time base');
            notXdf = EEG;
            notXdf.etc.alz.rawFile = fullfile(testCase.Folder, 'sub01.vhdr');
            testCase.verifySubstring(SyncTest.refusal(@() Heartbeats(notXdf, struct('Source', 'xdf'))), ...
                'not read from an XDF file');
        end

        % ---- route 'recording' -------------------------------------------
        function anotherRecordingsBeatsLandWhereTheyAre(testCase)
        %ANOTHERRECORDINGSBEATSLANDWHERETHEYARE  The ECG recorded at 1000 Hz on
        %   a clock 12.3 s off and drifting 40 ppm, three triggers missing there
        %   and two extra, found beside the EEG by its name.
            testCase.requireEeglab();
            [EEG, truth] = SyncTest.eeg();
            EEG.etc.alz.rawFile = fullfile(testCase.Folder, 'sub01.vhdr');
            [other, toEeg] = SyncTest.otherRecording(EEG, truth, 0.0001, [5 20 33], 2);
            pop_saveset(other, 'filename', 'sub01.set', 'filepath', testCase.Folder);
            for companion = {'sub01.eeg', 'sub01.vmrk'}         % the EEG's own BrainVision files
                fclose(fopen(fullfile(testCase.Folder, companion{1}), 'w'));
            end
            detected = toEeg(Rpeak.beats(double(other.data), other.srate, [], SyncTest.Detection));

            out = Heartbeats(EEG, struct('Source', 'recording'));
            beats = out.event(strcmp({out.event.type}, 'heartbeat'));
            testCase.assertNumElements(beats, numel(truth));
            testCase.verifyLessThan(max(abs([beats.latency] - detected)), 0.05, ...
                'Carried onto the EEG as the true clocks carry them, to 0.1 ms.');
            testCase.verifyLessThan(max(abs([beats.latency] - truth)), 1, ...
                'Within a sample (2 ms) of the true beats.');
            ibi = [beats.ibi];
            testCase.verifyEqual(ibi(1:end - 1), diff(detected) / SyncTest.Srate * 1000, 'AbsTol', 0.01, ...
                'The intervals, on the EEG''s clock.');
            sync = out.etc.alz.heartbeats.sync;
            testCase.verifyEqual(sync.route, 'recording');
            testCase.verifyEqual(sync.pairs, numel(EEG.event) - 3);
            testCase.verifyEqual(sync.unpairedHere, 3);
            testCase.verifyEqual(sync.unpairedThere, 2);
            testCase.verifyEqual(sync.driftPpm, 40, 'AbsTol', 1);
            testCase.verifyLessThan(sync.residualMaxMs, 0.2);
            testCase.verifyTrue(out.etc.alz.heartbeats.check.seen);
            testCase.verifyEqual(out.etc.alz.heartbeats.check.lagMs, 0, 'AbsTol', 1.5);
        end

        function recordingsThatDoNotLineUpAreRefused(testCase)
        %RECORDINGSTHATDONOTLINEUPAREREFUSED  Triggers timed up to 5 ms off at
        %   random; no code in common; a recording that is not there; a cut
        %   EEG; and triggers of one code every 2 s exactly, which fit an
        %   offset by a whole period nearly as well.
            testCase.requireEeglab();
            [EEG, truth] = SyncTest.eeg();
            EEG.etc.alz.rawFile = fullfile(testCase.Folder, 'sub01.vhdr');
            other = SyncTest.otherRecording(EEG, truth, 0.010, [], 0);
            pop_saveset(other, 'filename', 'sub01.set', 'filepath', testCase.Folder);
            refusal = @(varargin) SyncTest.refusal(@() Heartbeats(EEG, struct('Source', 'recording', varargin{:})));
            testCase.verifySubstring(refusal(), 'do not line up well enough');
            testCase.verifySubstring(refusal('SharedCodes', 'nothing'), 'no trigger code in common');
            testCase.verifySubstring(refusal('Recording', 'gone.set'), 'no recording gone.set');
            cut = EEG;
            cut.event(end + 1) = struct('type', 'boundary', 'latency', 1000.5, 'duration', 0);
            testCase.verifySubstring(SyncTest.refusal(@() Heartbeats(cut, struct('Source', 'recording'))), ...
                'has been cut');

            regular = EEG;
            regular.event = struct('type', '7', 'latency', num2cell(1 + (2:2:86) * SyncTest.Srate), 'duration', 0);
            regular.etc.alz.rawFile = fullfile(testCase.Folder, 'sub02.vhdr');
            other = SyncTest.otherRecording(regular, truth, 0.0001, [], 0);
            pop_saveset(other, 'filename', 'sub02.set', 'filepath', testCase.Folder);
            testCase.verifySubstring(SyncTest.refusal(@() Heartbeats(regular, struct('Source', 'recording'))), ...
                'too regular');
        end
    end

    methods (Access = private)
        function requireEeglab(testCase)
            if isempty(which('pop_saveset')) && ~isempty(which('eeglab'))
                evalc('eeglab(''nogui'');');
            end
            testCase.assumeNotEmpty(which('pop_saveset'), 'EEGLAB is not on the path.');
        end

        function requireXdf(testCase)
            testCase.requireEeglab();
            testCase.assumeNotEmpty(which('pop_loadxdf'), 'EEGLAB''s xdfimport plugin is not installed.');
            if isempty(which('load_xdf'))   % added by xdfimport's EEGLAB hook, which may not have run
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                    fullfile(fileparts(which('pop_loadxdf')), 'xdf')));
            end
        end

        function [EEG, truthLatency, detected, inGap] = xdfRecording(testCase, ecgRate, deviceDelay, gap)
        %XDFRECORDING  An XDF file with three streams, each on its own clock:
        %   the EEG (500 Hz, from recorder time 1010 s, the samples GAP lost),
        %   a chest strap's accelerometer (52 Hz) and its ECG (ECGRATE, from
        %   1005 s, 10 s longer at each end, its time stamps DEVICEDELAY
        %   late); and the EEG read from it as Alakazam reads it.
        %   TRUTHLATENCY: every beat's place in the EEG, in samples (with no
        %   gap). DETECTED: where the true clocks put the beats the detector
        %   finds in the ECG as stamped, NaN outside the EEG and in its gap;
        %   INGAP: which of them fell in the gap.
            if nargin < 4
                gap = [];
            end
            srate = SyncTest.Srate;
            stream = RandStream('mt19937ar', 'Seed', 21);
            eegStart = 1010;
            ecgStart = 1005;
            beatTimes = SyncTest.beats(ecgStart + 0.6, ecgStart + SyncTest.Seconds + 19, stream);
            n = SyncTest.Seconds * srate;
            eegTrue = eegStart + (0:n - 1) / srate;
            ecgTrue = ecgStart + (0:round((SyncTest.Seconds + 20) * ecgRate) - 1) / ecgRate;
            ecg = SyncTest.ecgAt(ecgTrue, beatTimes) + 10 * randn(stream, 1, numel(ecgTrue));
            field = SyncTest.ecgAt(eegTrue, beatTimes);
            eegData = [0.05; -0.03; 0.04; 0.02] .* field + 2 * randn(stream, 4, n);
            kept = setdiff(1:n, gap);
            accTrue = ecgStart + (0:round((SyncTest.Seconds + 20) * 52) - 1) / 52;

            eegOffset = 1.25;
            ecgOffset = -3.7;
            jitter = @(m) 2e-4 * (rand(stream, 1, m) - 0.5);
            offsets = @(value) [(980:5:1120)', value * ones(29, 1)];
            streams = struct( ...
                'name', {'Amp EEG', 'Chest ACC', 'Chest ECG'}, 'type', {'EEG', 'Accelerometer', 'ECG'}, ...
                'labels', {{'Fz', 'Cz', 'Pz', 'Oz'}, {'X', 'Y', 'Z'}, {'ECG'}}, ...
                'rate', {srate, 52, ecgRate}, ...
                'times', {eegTrue(kept) - eegOffset + jitter(numel(kept)), ...
                          accTrue - ecgOffset + deviceDelay + jitter(numel(accTrue)), ...
                          ecgTrue - ecgOffset + deviceDelay + jitter(numel(ecgTrue))}, ...
                'data', {eegData(:, kept), randn(stream, 3, numel(accTrue)), ecg}, ...
                'offsets', {offsets(eegOffset), offsets(ecgOffset), offsets(ecgOffset)});
            file = fullfile(testCase.Folder, 'sub01.xdf');
            SyncTest.writeXdf(file, streams);

            [~, EEG] = evalc('pop_loadxdf(file, ''streamtype'', ''EEG'')');
            EEG = eeg_checkset(EEG);
            EEG = recordRawFile(EEG, file);
            EEG.DataType = 'TIMEDOMAIN';
            EEG.DataFormat = 'CONTINUOUS';
            truthLatency = 1 + (beatTimes - eegStart) * srate;
            found = Rpeak.beats(double(single(ecg)), ecgTrue + deviceDelay, [], SyncTest.Detection);
            detected = interp1(eegTrue(kept), 1:numel(kept), found, 'linear', NaN);
            inGap = false(size(found));
            if ~isempty(gap)
                inGap = found > eegTrue(gap(1) - 1) & found < eegTrue(gap(end) + 1);
                detected(inGap) = NaN;
            end
        end
    end

    methods (Static)
        function text = refusal(fn)
        %REFUSAL  What Heartbeats said when it refused; '' when it did not.
            text = '';
            try
                fn();
            catch err
                if ~strcmp(err.identifier, 'Alakazam:Heartbeats')
                    rethrow(err);
                end
                text = err.message;
            end
        end

        function ms = lateBy(EEG, truthLatency)
        %LATEBY  How late each heartbeat event lies after the true beat
        %   nearest it, in ms.
            latency = [EEG.event(strcmp({EEG.event.type}, 'heartbeat')).latency];
            [~, nearest] = min(abs(latency - truthLatency.'), [], 1);
            ms = (latency - truthLatency(nearest)) / EEG.srate * 1000;
        end

        function times = beats(from, to, stream)
        %BEATS  Heartbeat times from FROM to TO: about 75 a minute, varying.
            times = from;
            while times(end) < to
                times(end + 1) = times(end) + 0.8 + 0.06 * sin(2 * pi * 0.1 * times(end)) ...
                    + 0.02 * randn(stream); %#ok<AGROW>
            end
            times = times(1:end - 1);
        end

        function ecg = ecgAt(t, beatTimes)
        %ECGAT  An ECG at times T (s) with R-peaks at BEATTIMES: P, Q, R, S, T
        %   waves, in microvolts.
            ecg = zeros(size(t));
            for b = beatTimes(beatTimes > t(1) - 1 & beatTimes < t(end) + 1)
                near = abs(t - b) < 0.6;
                tt = t(near) - b;
                ecg(near) = ecg(near) + 100 * exp(-0.5 * ((tt + 0.16) / 0.025) .^ 2) ...
                    - 120 * exp(-0.5 * ((tt + 0.025) / 0.006) .^ 2) ...
                    + 1000 * exp(-0.5 * (tt / 0.008) .^ 2) ...
                    - 200 * exp(-0.5 * ((tt - 0.025) / 0.006) .^ 2) ...
                    + 150 * exp(-0.5 * ((tt - 0.26) / 0.04) .^ 2);
            end
        end

        function [EEG, truth] = eeg()
        %EEG  A continuous EEG at 500 Hz, unfiltered: four scalp channels,
        %   each at its own DC level, carrying the heart's field; and 40
        %   triggers of codes 1 to 3 at irregular times. TRUTH: every beat's
        %   latency in samples.
            srate = SyncTest.Srate;
            stream = RandStream('mt19937ar', 'Seed', 11);
            n = SyncTest.Seconds * srate;
            t = (0:n - 1) / srate;
            beatTimes = SyncTest.beats(0.7, SyncTest.Seconds - 1, stream);
            field = SyncTest.ecgAt(t, beatTimes);
            EEG = struct();
            EEG.srate = srate;
            EEG.nbchan = 4;
            EEG.pnts = n;
            EEG.trials = 1;
            EEG.xmin = 0;
            EEG.xmax = (n - 1) / srate;
            EEG.times = t;
            offsets = [300; -200; 150; 80];                     % unfiltered: each channel at its own level
            EEG.data = single(offsets + [0.05; -0.03; 0.04; 0.02] .* field + 2 * randn(stream, 4, n));
            EEG.chanlocs = struct('labels', {'Fz', 'Cz', 'Pz', 'Oz'}, 'type', {'EEG', 'EEG', 'EEG', 'EEG'});
            EEG.DataFormat = 'CONTINUOUS';
            EEG.DataType = 'TIMEDOMAIN';
            EEG.etc = struct();
            EEG.event = struct('type', arrayfun(@(k) sprintf('%d', mod(k, 3) + 1), 1:40, 'UniformOutput', false), ...
                'latency', num2cell(round(sort(1.5 + (SyncTest.Seconds - 3) * rand(stream, 1, 40)) * srate) + 1), ...
                'duration', 0);
            EEG.setname = 'made up';
            truth = 1 + beatTimes * srate;
        end

        function [other, toEeg] = otherRecording(EEG, truth, triggerJitter, missing, extra)
        %OTHERRECORDING  The ECG recorded by another device at 1000 Hz, on a
        %   clock 12.3 s off the EEG's and drifting 40 ppm from it, with the
        %   EEG's triggers (each timed to within TRIGGERJITTER s, those in
        %   MISSING left out, EXTRA more of code 1 at random times), as an
        %   EEGLAB dataset. TOEEG carries a time in it (s from its first
        %   sample) to its latency in the EEG, as the true clocks do.
            stream = RandStream('mt19937ar', 'Seed', 5);
            rate = 1000;
            a = 12.3;
            b = 1 + 40e-6;
            beatTimes = (((truth - 1) / EEG.srate) - a) / b;          % on the other clock
            tEeg = ([EEG.event.latency] - 1) / EEG.srate;
            tThere = (tEeg - a) / b + triggerJitter * (rand(stream, 1, numel(tEeg)) - 0.5);
            codes = {EEG.event.type};
            tThere(missing) = [];
            codes(missing) = [];
            tThere = [tThere, -a + 5 + 20 * rand(stream, 1, extra)];
            codes = [codes, repmat({'1'}, 1, extra)];
            start = -a - 5;                                            % its first sample's time
            n = round((SyncTest.Seconds + 10) * rate);
            t = start + (0:n - 1) / rate;
            other = eeg_emptyset();
            other.data = single(SyncTest.ecgAt(t, beatTimes) + 10 * randn(stream, 1, n));
            other.srate = rate;
            other.nbchan = 1;
            other.pnts = n;
            other.trials = 1;
            other.xmin = 0;
            other.xmax = (n - 1) / rate;
            other.chanlocs = struct('labels', 'ECG', 'type', 'ECG');
            other.event = struct('type', codes, 'latency', num2cell(1 + (tThere - start) * rate), ...
                'duration', 0);
            other.setname = 'other device';
            toEeg = @(time) 1 + (a + b * (start + time)) * EEG.srate;
            other = eeg_checkset(other, 'eventconsistency');
        end

        function writeXdf(file, streams)
        %WRITEXDF  A minimal XDF 1.0 file: the file header, each stream's
        %   header, its samples (float32, each with its time stamp), its clock
        %   offsets [collection time, offset] and its footer, which lists the
        %   clock offsets again, as LabRecorder writes it.
            fid = fopen(file, 'w', 'l');
            closeFile = onCleanup(@() fclose(fid));
            fwrite(fid, 'XDF:', 'char');
            chunk(fid, 1, uint8('<?xml version="1.0"?><info><version>1.0</version></info>'));
            for s = 1:numel(streams)
                S = streams(s);
                channels = strjoin(cellfun(@(l) sprintf(['<channel><label>%s</label><type>%s</type>' ...
                    '<unit>microvolts</unit></channel>'], l, S.type), S.labels, 'UniformOutput', false), '');
                xml = sprintf(['<?xml version="1.0"?><info><name>%s</name><type>%s</type>' ...
                    '<channel_count>%d</channel_count><nominal_srate>%g</nominal_srate>' ...
                    '<channel_format>float32</channel_format><source_id>%s</source_id>' ...
                    '<version>1.1</version><created_at>0</created_at><uid>u%d</uid>' ...
                    '<session_id>default</session_id><hostname>test</hostname>' ...
                    '<desc><channels>%s</channels></desc></info>'], S.name, S.type, numel(S.labels), ...
                    S.rate, S.name, s, channels);
                chunk(fid, 2, [typecast(uint32(s), 'uint8'), uint8(xml)]);
            end
            for s = 1:numel(streams)
                S = streams(s);
                n = numel(S.times);
                channelCount = size(S.data, 1);
                for first = 1:1000:n
                    k = first:min(n, first + 999);
                    perSample = [repmat(uint8(8), 1, numel(k)); ...
                        reshape(typecast(double(S.times(k)), 'uint8'), 8, []); ...
                        reshape(typecast(reshape(single(S.data(:, k)), 1, []), 'uint8'), 4 * channelCount, [])];
                    chunk(fid, 3, [typecast(uint32(s), 'uint8'), uint8(4), ...
                        typecast(uint32(numel(k)), 'uint8'), reshape(perSample, 1, [])]);
                end
                for c = 1:size(S.offsets, 1)
                    chunk(fid, 4, [typecast(uint32(s), 'uint8'), typecast(double(S.offsets(c, 1)), 'uint8'), ...
                        typecast(double(S.offsets(c, 2)), 'uint8')]);
                end
            end
            for s = 1:numel(streams)
                S = streams(s);
                offsets = strjoin(arrayfun(@(c) sprintf('<offset><time>%.10f</time><value>%.10f</value></offset>', ...
                    S.offsets(c, 1), S.offsets(c, 2)), 1:size(S.offsets, 1), 'UniformOutput', false), '');
                xml = sprintf(['<?xml version="1.0"?><info><first_timestamp>%.10f</first_timestamp>' ...
                    '<last_timestamp>%.10f</last_timestamp><sample_count>%d</sample_count>' ...
                    '<clock_offsets>%s</clock_offsets></info>'], S.times(1), S.times(end), numel(S.times), offsets);
                chunk(fid, 6, [typecast(uint32(s), 'uint8'), uint8(xml)]);
            end
        end
    end
end

% ======================================================================= %
function chunk(fid, tag, content)
%CHUNK  One XDF chunk: its length in 8 bytes, its tag, its content.
    fwrite(fid, uint8(8), 'uint8');
    fwrite(fid, uint64(numel(content) + 2), 'uint64');
    fwrite(fid, uint16(tag), 'uint16');
    fwrite(fid, content, 'uint8');
end
