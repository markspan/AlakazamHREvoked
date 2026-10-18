classdef HeartbeatsTest < matlab.unittest.TestCase
%HEARTBEATSTEST  Heartbeats on a made-up recording whose heartbeats are known:
%   three EEG channels of noise and an ECG of 140 beats at a varying rate,
%   with baseline wander, P and T waves and noise, and a stimulus every two
%   seconds, alternately 1 and 2.
%
%   The beats are found where they are, the same when the ECG is inverted,
%   separately on either side of a cut or a rejected stretch, and not twice
%   on a second run; the data are untouched; DefineBins chooses them by the
%   stimulus before them; the ECG channel is found by its type or label;
%   and what cannot work is refused. That the detector is spectHR's is
%   SpectHRParityTest's.
%
%   Run with: runtests('tests/HeartbeatsTest.m').
%
%   See also HEARTBEATS, SPECTHRPARITYTEST.

    properties (Constant)
        Srate = 500
        Seconds = 120
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
                     fullfile(root, 'src', 'Transformations', 'DefineBins'), ...
                     fullfile(root, 'src', 'Support'), fullfile(root, 'src', 'Dialogs'), ...
                     fullfile(plugin, 'Heartbeats')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(p{1}));
            end
        end
    end

    methods (Test)
        function theBeatsAreFoundWhereTheyAre(testCase)
            [EEG, truth] = HeartbeatsTest.recording();
            out = Heartbeats(EEG, struct('Channel', 'ECG'));

            beats = HeartbeatsTest.beats(out);
            testCase.assertNumElements(beats, numel(truth), 'One event per heartbeat.');
            testCase.verifyLessThanOrEqual(max(abs([beats.latency] - truth)), 0.5, ...
                'Each R-peak within half a sample of where it is.');
            ibi = [beats.ibi];
            testCase.verifyEqual(ibi(1:end - 1), diff(truth) / HeartbeatsTest.Srate * 1000, 'AbsTol', 2, ...
                'The interval to the next beat, in ms.');
            testCase.verifyTrue(isnan(ibi(end)), 'The last beat has no next one.');
            testCase.verifyEqual([beats.hr], 60000 ./ ibi);
            testCase.verifyEqual(out.data, EEG.data, 'The data are not changed.');

            stimuli = out.event(~strcmp({out.event.type}, 'heartbeat'));
            testCase.verifyEqual({stimuli.type}, {EEG.event.type}, 'The other events stay.');
            testCase.verifyEqual([stimuli.latency], [EEG.event.latency]);
            testCase.verifyTrue(issorted([out.event.latency]), 'The events stay in time order.');

            record = out.etc.alz.heartbeats;
            testCase.verifyEqual(record.channel, 'ECG');
            testCase.verifyEqual(record.polarity, 'normal');
            testCase.verifyEqual(record.polarityBy, 'skewness');
            testCase.verifyEqual(record.beats, numel(truth));
            testCase.verifyEqual(record.stretches, 1);
            testCase.verifyEqual(record.meanHr, mean(60000 ./ ibi(1:end - 1)), 'RelTol', 1e-12);
        end

        function anInvertedEcgIsTurnedUprightToFindThem(testCase)
            [EEG, truth] = HeartbeatsTest.recording();
            ecg = strcmp({EEG.chanlocs.labels}, 'ECG');
            EEG.data(ecg, :) = -EEG.data(ecg, :);
            out = Heartbeats(EEG, struct('Channel', 'ECG'));

            beats = HeartbeatsTest.beats(out);
            testCase.assertNumElements(beats, numel(truth));
            testCase.verifyLessThanOrEqual(max(abs([beats.latency] - truth)), 0.5);
            testCase.verifyEqual(out.etc.alz.heartbeats.polarity, 'inverted');
            testCase.verifyEqual(out.data, EEG.data, 'Turned upright to search, not in the data.');

            told = Heartbeats(EEG, struct('Channel', 'ECG', 'Polarity', 'Inverted'));
            testCase.verifyEqual([HeartbeatsTest.beats(told).latency], [beats.latency]);
            testCase.verifyEqual(told.etc.alz.heartbeats.polarityBy, 'the user');
        end

        function cutsAndRejectedStretchesAreSearchedApart(testCase)
        %CUTSANDREJECTEDSTRETCHESARESEARCHEDAPART  A cut near 60 s and a
        %   rejected (NaN) stretch near 90 to 92 s, each placed midway between
        %   two beats: three stretches, no beat in the rejected one, and the
        %   last beat before each gap without an interval, since the next
        %   beat is on the other side.
            [EEG, truth] = HeartbeatsTest.recording();
            srate = HeartbeatsTest.Srate;
            midway = @(seconds) round(mean(truth(find(truth < seconds * srate, 1, 'last') + [0 1])));
            cut = midway(60);
            EEG.event(end + 1) = struct('type', 'boundary', 'latency', cut + 0.5, 'duration', 0);
            [~, order] = sort([EEG.event.latency]);
            EEG.event = EEG.event(order);
            rejected = midway(90):midway(92);
            EEG.data(:, rejected) = NaN;
            out = Heartbeats(EEG, struct('Channel', 'ECG'));

            beats = HeartbeatsTest.beats(out);
            latency = [beats.latency];
            testCase.verifyEqual(out.etc.alz.heartbeats.stretches, 3);
            testCase.verifyFalse(any(latency >= rejected(1) & latency <= rejected(end)), ...
                'No beat in a rejected stretch.');
            kept = truth(~(truth >= rejected(1) & truth <= rejected(end)));
            testCase.assertNumElements(latency, numel(kept));
            testCase.verifyLessThanOrEqual(max(abs(latency - kept)), 0.5);
            ibi = [beats.ibi];
            beforeCut = find(latency < cut, 1, 'last');
            beforeGap = find(latency < rejected(1), 1, 'last');
            testCase.verifyTrue(isnan(ibi(beforeCut)), 'No interval across a cut.');
            testCase.verifyTrue(isnan(ibi(beforeGap)), 'No interval across a rejected stretch.');
            testCase.verifyEqual(nnz(isnan(ibi)), 3, 'Only the last beat of each stretch.');
        end

        function aSecondRunReplacesTheBeats(testCase)
            [EEG, truth] = HeartbeatsTest.recording();
            once = Heartbeats(EEG, struct('Channel', 'ECG'));
            twice = Heartbeats(once, struct('Channel', 'ECG'));

            testCase.verifyNumElements(HeartbeatsTest.beats(twice), numel(truth));
            testCase.verifyEqual(twice.etc.alz.heartbeats.replaced, numel(truth));
            testCase.verifyEqual(twice.event, once.event);
        end

        function definesBinsChooseThemByWhatCameBefore(testCase)
        %DEFINESBINSCHOOSETHEMBYWHATCAMEBEFORE  The heartbeats in the second
        %   after a 1, and in the second after a 2, as the manual's example.
            EEG = HeartbeatsTest.recording();
            out = Heartbeats(EEG, struct('Channel', 'ECG'));
            script = sprintf(['bin 1 "After 1" "heartbeat" any(1) within [-1000,0] ms\n' ...
                'bin 2 "After 2" "heartbeat" any(2) within [-1000,0] ms']);
            tagged = DefineBins(out, struct('script', script));

            beat = find(strcmp({tagged.event.type}, 'heartbeat'));
            srate = HeartbeatsTest.Srate;
            for bin = 1:2
                code = sprintf('%d', bin);
                stimuli = [tagged.event(strcmp({tagged.event.type}, code)).latency];
                expected = beat(arrayfun(@(k) any(tagged.event(k).latency - stimuli >= 0 & ...
                    tagged.event(k).latency - stimuli <= srate), beat));
                testCase.verifyEqual(tagged.bindesc(bin).events, expected, sprintf('Bin %d.', bin));
                testCase.verifyGreaterThan(numel(expected), 20);
            end
        end

        function theEcgChannelIsFoundByItsTypeOrLabel(testCase)
            [EEG, truth] = HeartbeatsTest.recording();
            out = Heartbeats(EEG, struct());
            testCase.verifyEqual(out.etc.alz.heartbeats.channel, 'ECG', 'By its label.');

            EEG.chanlocs(4).labels = 'BIP3';
            [EEG.chanlocs.type] = deal('EEG');
            EEG.chanlocs(4).type = 'ECG';
            out = Heartbeats(EEG, struct());
            testCase.verifyEqual(out.etc.alz.heartbeats.channel, 'BIP3', 'By its type.');
            testCase.verifyNumElements(HeartbeatsTest.beats(out), numel(truth));
        end

        function whatCannotWorkIsRefused(testCase)
            EEG = HeartbeatsTest.recording();
            testCase.verifyError(@() Heartbeats(EEG, struct('Channel', 'ECG2')), 'Alakazam:Heartbeats');
            noEcg = EEG;
            [noEcg.chanlocs.labels] = deal('Fz', 'Cz', 'Pz', 'Oz');
            testCase.verifyError(@() Heartbeats(noEcg, struct()), 'Alakazam:Heartbeats');
            flat = EEG;
            flat.data(4, :) = 0;
            testCase.verifyError(@() Heartbeats(flat, struct('Channel', 'ECG')), 'Alakazam:Heartbeats');
            epoched = EEG;
            epoched.DataFormat = 'EPOCHED';
            testCase.verifyError(@() Heartbeats(epoched, struct('Channel', 'ECG')), 'Alakazam:Heartbeats');
            testCase.verifyError(@() Heartbeats(EEG, struct('Channel', 'ECG', 'HighpassHz', 300)), ...
                'Alakazam:Heartbeats');
        end

        function theDialogOffersTheEcgChannelAndSpectHRsSettings(testCase)
        %THEDIALOGOFFERSTHEECGCHANNELANDSPECTHRSSETTINGS  Opened on the
        %   recording and OK pressed: the ECG from this recording's channel,
        %   the other sources' fields greyed out, spectHR's defaults and the
        %   check on. A timer presses OK, or Cancel when OK is refused.
            try
                probe = uifigure('Visible', 'off');
                delete(probe);
            catch err
                testCase.assumeFail(['A uifigure could not be created here: ' err.message]);
            end
            EEG = HeartbeatsTest.recording();
            state = struct('blankText', [], 'misfit', []);
            pressTimer = timer('StartDelay', 8, 'TimerFcn', @(~, ~) press());
            cleanup = onCleanup(@() delete(pressTimer));
            start(pressTimer);
            [out, options] = Heartbeats(EEG);
            clear cleanup;

            testCase.assertNotEmpty(options, 'The dialog was cancelled.');
            testCase.verifyEqual(options.Source, 'channel');
            testCase.verifyEqual(options.Channel, 'ECG');
            testCase.verifyEqual(options.Polarity, 'Automatic');
            testCase.verifyEqual(options.MinIntervalMs, 300);
            testCase.verifyEqual(options.HighpassHz, 0.5);
            testCase.verifyEqual(options.EventType, 'heartbeat');
            testCase.verifyTrue(options.Check);
            testCase.verifyFalse(options.CorrectLag);
            testCase.verifyNumElements(HeartbeatsTest.beats(out), 140);
            off = matlab.lang.OnOffSwitchState('off');
            testCase.verifyEqual(state.blankText, repmat(off, 1, 5), ...
                'Stream, its channel, recording, its channel, shared codes: greyed out.');
            testCase.verifyEqual(state.misfit, off, 'The largest misfit: greyed out.');

            function press()
                f = findall(groot, 'Type', 'figure', 'Name', 'Heartbeats');
                if isempty(f)
                    return;
                end
                texts = findall(f(1), 'Type', 'uieditfield');
                state.blankText = [texts(strcmp(get(texts, {'Value'}), '')).Enable];
                numbers = findall(f(1), 'Type', 'uinumericeditfield');
                state.misfit = [numbers([numbers.Value] == 2).Enable];
                ok = findall(f(1), 'Type', 'uibutton', 'Text', 'OK');
                if ~isempty(ok)
                    ok(1).ButtonPushedFcn(ok(1), []);
                end
                if isvalid(f(1))
                    cancel = findall(f(1), 'Type', 'uibutton', 'Text', 'Cancel');
                    cancel(1).ButtonPushedFcn(cancel(1), []);
                end
            end
        end
    end

    methods (Static)
        function [EEG, truth] = recording()
        %RECORDING  The made-up recording, and its R-peaks' latencies in
        %   samples (1-based, fractional).
            srate = HeartbeatsTest.Srate;
            n = HeartbeatsTest.Seconds * srate;
            t = (0:n - 1) / srate;
            stream = RandStream('mt19937ar', 'Seed', 11);

            beatTimes = 0.7;
            while beatTimes(end) < HeartbeatsTest.Seconds - 1.2
                rr = 0.8 + 0.06 * sin(2 * pi * 0.1 * beatTimes(end)) + 0.02 * randn(stream);
                beatTimes(end + 1) = beatTimes(end) + rr; %#ok<AGROW>
            end
            beatTimes = beatTimes(1:140);

            wave = @(centre, width, height) height * exp(-0.5 * ((t - centre) / width) .^ 2);
            ecg = 300 * sin(2 * pi * 0.2 * t) + 10 * randn(stream, 1, n);
            for b = beatTimes
                ecg = ecg + wave(b - 0.16, 0.025, 100) + wave(b - 0.025, 0.006, -120) ...
                    + wave(b, 0.008, 1000) + wave(b + 0.025, 0.006, -200) + wave(b + 0.26, 0.04, 250);
            end

            EEG = struct();
            EEG.srate = srate;
            EEG.nbchan = 4;
            EEG.pnts = n;
            EEG.trials = 1;
            EEG.xmin = 0;
            EEG.xmax = (n - 1) / srate;
            EEG.times = t;
            EEG.data = single([20 * randn(stream, 3, n); ecg]);
            EEG.chanlocs = struct('labels', {'Fz', 'Cz', 'Pz', 'ECG'}, 'type', {'', '', '', ''});
            EEG.DataFormat = 'CONTINUOUS';
            EEG.DataType = 'TIMEDOMAIN';
            EEG.etc = struct();
            onsets = 1:2:HeartbeatsTest.Seconds - 2;
            codes = repmat({'1', '2'}, 1, ceil(numel(onsets) / 2));
            EEG.event = struct('type', codes(1:numel(onsets)), ...
                'latency', num2cell(round(onsets * srate + 100 * rand(stream, 1, numel(onsets))) + 1), ...
                'duration', 0);
            truth = beatTimes * srate + 1;
        end

        function beats = beats(EEG)
            beats = EEG.event(strcmp({EEG.event.type}, 'heartbeat'));
        end
    end
end
