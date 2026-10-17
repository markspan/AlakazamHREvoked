classdef SpectHRParityTest < matlab.unittest.TestCase
%SPECTHRPARITYTEST  The R-peaks are spectHR's: every step of the MATLAB port
%   (the +Rpeak package) against spectHR's own results on the same signals.
%
%   fixtures/spectHRReference.mat was made by spectHR's own code
%   (fixtures/makeSpectHRReference.py, run with spectHR's Python): five
%   minutes of the ECG of spectHR's VU-AMS example recording, its polarity,
%   its 0.5 Hz high-pass and its R-peaks; the same at exactly 500 Hz, as an
%   EEG recording stores its rate; two made-up signals that reach the
%   polarity tiebreaker; and scipy's peak finder on plateaus and close
%   peaks. numpy's sums must come out to the last bit, and so must the peak
%   finder, the refractory distance in samples (299 at this recording's
%   rate, 999.9999998 Hz) and, on the unfiltered ECG, the threshold and the
%   beats' times; after the filter, the R-peaks must be the same beats on
%   the same samples, their times within 1e-10 s.
%
%   THE FILTERS TO A MILLIONTH. numpy multiplies complex numbers with fused
%   multiply-adds (in np.prod, and in the BLAS that np.poly's convolutions
%   call), which MATLAB cannot reproduce; whether numpy does depends on the
%   processor, so spectHR's own coefficients differ in the last bit from one
%   machine to another. The coefficients here differ from scipy's in their
%   last one or two bits, and a 0.5 Hz high-pass at 1000 Hz, its poles a few
%   thousandths from 1, magnifies that to a few parts in ten million of the
%   signal (7e-7 on a signal of +-3.1). Smith's division, which numpy also
%   uses, is reproduced exactly. The R-peaks do not move.
%
%   Run with: runtests('tests/SpectHRParityTest.m').
%
%   See also RPEAK.DETECT, RPEAK.POLARITY, RPEAK.HIGHPASS, RPEAK.NPSUM.

    properties (Constant)
        Plugin = fileparts(fileparts(mfilename('fullpath')))
    end

    properties
        Ref
    end

    methods (TestClassSetup)
        function load(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(SpectHRParityTest.Plugin, 'Heartbeats')));
            testCase.Ref = builtin('load', fullfile(SpectHRParityTest.Plugin, 'tests', 'fixtures', ...
                'spectHRReference.mat'));
        end
    end

    methods (Test)
        function numpysSumsToTheLastBit(testCase)
            ref = testCase.Ref;
            times = (0:numel(ref.ecg) - 1) / ref.srate;
            testCase.verifyEqual(Rpeak.npSum(ref.ecg), ref.numpySum);
            testCase.verifyEqual(Rpeak.npMean(diff(times)), ref.numpyMeanStep);
            testCase.verifyEqual(Rpeak.npStd(ref.ecg), ref.numpyStd);
            testCase.verifyEqual(Rpeak.npMedian(ref.ecg), ref.numpyMedian);
        end

        function thePolarityIsSpectHRs(testCase)
            ref = testCase.Ref;
            times = (0:numel(ref.ecg) - 1) / ref.srate;
            [polarity, info] = Rpeak.polarity(times, ref.ecg);
            testCase.verifyEqual(polarity, ref.polarity);
            testCase.verifyEqual(info.skewness, ref.skewness, 'RelTol', 1e-6, ...
                'Band-passed through a filter a bit or two off scipy''s (see the class help).');
            [polarity, info] = Rpeak.polarity(times, -ref.ecg);
            testCase.verifyEqual(polarity, ref.invertedPolarity);
            testCase.verifyEqual(info.skewness, ref.invertedSkewness, 'RelTol', 1e-6);
        end

        function theTiebreakerIsSpectHRs(testCase)
        %THETIEBREAKERISSPECTHRS  Where the skewness is too small to trust, the
        %   peak prominences decide, once for an inverted and once for an
        %   upright made-up signal.
            ref = testCase.Ref;
            for name = {'tie', 'tieUp'}
                [polarity, info] = Rpeak.polarity(ref.tieTimes, ref.(name{1}));
                testCase.verifyEqual(info.decidedBy, 'peak_count', name{1});
                testCase.verifyEqual(polarity, ref.([name{1} 'Polarity']), name{1});
                testCase.verifyEqual(info.skewness, ref.([name{1} 'Skewness']), 'RelTol', 1e-9, name{1});
                testCase.verifyEqual(info.peakScore, ref.([name{1} 'Score']), 'RelTol', 1e-9, name{1});
            end
        end

        function theHighpassIsSpectHRs(testCase)
        %THEHIGHPASSISSPECTHRS  At its ends, where scipy's filtfilt differs
        %   from MATLAB's, and in its middle.
            ref = testCase.Ref;
            times = (0:numel(ref.ecg) - 1) / ref.srate;
            filtered = Rpeak.highpass(ref.ecg, times, 0.5);
            testCase.verifyEqual(filtered(ref.highpassedAt), ref.highpassed, ...
                'AbsTol', 1e-6 * max(abs(ref.ecg)), 'See the class help: a millionth.');
        end

        function theRpeaksAreSpectHRs(testCase)
            ref = testCase.Ref;
            times = (0:numel(ref.ecg) - 1) / ref.srate;
            [peakTimes, ~, info] = Rpeak.detect(times, Rpeak.highpass(ref.ecg, times, 0.5), 300);
            testCase.verifyEqual(info.minDistanceSamples, double(ref.minDistanceSamples), ...
                'The refractory distance in samples, as spectHR truncates it.');
            testCase.assertNumElements(peakTimes, numel(ref.peakTimes));
            testCase.verifyEqual(peakTimes, ref.peakTimes, 'AbsTol', 1e-10);
        end

        function theDetectorAloneToTheLastBit(testCase)
        %THEDETECTORALONETOTHELASTBIT  On the unfiltered ECG, where nothing
        %   but numpy's sums and plain arithmetic is involved: the threshold
        %   (the median plus 1.5 population standard deviations) and every
        %   beat's time exactly as spectHR's.
            ref = testCase.Ref;
            times = (0:numel(ref.ecg) - 1) / ref.srate;
            [peakTimes, ~, info] = Rpeak.detect(times, ref.ecg, 300);
            testCase.verifyEqual(info.threshold, ref.rawThreshold);
            testCase.assertNumElements(peakTimes, numel(ref.rawPeakTimes));
            testCase.verifyEqual(peakTimes, ref.rawPeakTimes);
        end

        function theRpeaksAtFiveHundredHertzAreSpectHRs(testCase)
            ref = testCase.Ref;
            ecg = ref.ecg(1:2:120000);
            times = (0:numel(ecg) - 1) / 500;
            polarity = Rpeak.polarity(times, ecg);
            testCase.verifyEqual(polarity, ref.polarity500);
            [peakTimes, ~, info] = Rpeak.detect(times, Rpeak.highpass(ecg, times, 0.5), 300);
            testCase.verifyEqual(info.minDistanceSamples, 150);
            testCase.assertNumElements(peakTimes, numel(ref.peakTimes500));
            testCase.verifyEqual(peakTimes, ref.peakTimes500, 'AbsTol', 1e-10);
        end

        function scipysPeakFinder(testCase)
        %SCIPYSPEAKFINDER  Plateaus of two and three samples, peaks closer than
        %   the distance, by height and distance, and by distance and
        %   prominence, with the prominences.
            ref = testCase.Ref;
            x = ref.findPeaksSignal;
            testCase.verifyEqual(Rpeak.findPeaks(x, 'Height', 0.5, 'Distance', 7), ...
                double(ref.findPeaksHeightDistance));
            [peaks, prominences] = Rpeak.findPeaks(x, 'Distance', 5, 'Prominence', 1.0);
            testCase.verifyEqual(peaks, double(ref.findPeaksDistanceProminence));
            testCase.verifyEqual(prominences, ref.prominences, 'AbsTol', 1e-12);
        end
    end
end
