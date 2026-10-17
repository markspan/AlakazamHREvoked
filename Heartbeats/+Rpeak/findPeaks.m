function [peaks, prominences] = findPeaks(x, varargin)
%FINDPEAKS  scipy.signal.find_peaks, with the conditions spectHR uses.
%   [PEAKS, PROMINENCES] = Rpeak.findPeaks(X, Name, Value) are the indices
%   (1-based) of the local maxima of X that meet the conditions, in the order
%   scipy applies them:
%     'Height'      a peak's value at least this
%     'Distance'    of peaks closer than this many samples, the higher one
%                   stays: the highest first, each removing its lower
%                   neighbours within the distance
%     'Prominence'  a peak's prominence at least this (computed after the
%                   distance, on the peaks it kept), the height above the
%                   higher of its two lowest points before the signal rises
%                   above it again, as scipy.signal.peak_prominences
%   PROMINENCES are the kept peaks' prominences ([] without 'Prominence').
%
%   A LOCAL MAXIMUM is a sample higher than the one before it and than the
%   next one that differs from it; a flat top counts once, at its middle (the
%   left one of two middles). The first and the last sample never are.
%
%   EQUAL HEIGHTS. Of two peaks of exactly equal height within the distance,
%   scipy keeps the one its (unstable) sort puts last; MATLAB's sort is
%   stable, so here it is the later one. Real ECG has no such ties.
%
%   See also RPEAK.DETECT, RPEAK.POLARITY.
    p = inputParser();
    p.addParameter('Height', []);
    p.addParameter('Distance', []);
    p.addParameter('Prominence', []);
    p.parse(varargin{:});
    o = p.Results;
    x = double(reshape(x, 1, []));
    prominences = [];

    peaks = localMaxima(x);
    if ~isempty(o.Height)
        peaks = peaks(o.Height <= x(peaks));
    end
    if ~isempty(o.Distance)
        peaks = peaks(byDistance(peaks, x(peaks), o.Distance));
    end
    if ~isempty(o.Prominence)
        prominences = peakProminences(x, peaks);
        keep = o.Prominence <= prominences;
        peaks = peaks(keep);
        prominences = prominences(keep);
    end
end

% ======================================================================= %
function peaks = localMaxima(x)
%LOCALMAXIMA  scipy's _local_maxima_1d, the midpoints of every local
%   maximum and flat top. Candidates first (a rise into the sample), then
%   scipy's walk over each candidate's flat top.
    n = numel(x);
    peaks = zeros(1, 0);
    if n < 3
        return;
    end
    rises = find(x(1:n - 2) < x(2:n - 1)) + 1;   % i with x(i-1) < x(i), i in 2..n-1
    count = 0;
    peaks = zeros(1, numel(rises));
    skipTo = 0;   % scipy moves past a flat top it has walked
    for i = rises
        if i <= skipTo
            continue;
        end
        ahead = i + 1;
        while ahead < n && x(ahead) == x(i)
            ahead = ahead + 1;
        end
        if x(ahead) < x(i)
            count = count + 1;
            peaks(count) = floor((i - 1 + ahead - 2) / 2) + 1;   % 0-based (left + right) // 2
            skipTo = ahead;
        end
    end
    peaks = peaks(1:count);
end

function keep = byDistance(peaks, priority, distance)
%BYDISTANCE  scipy's _select_by_peak_distance.
    n = numel(peaks);
    distance = ceil(distance);
    keep = true(1, n);
    [~, order] = sort(priority);
    for i = n:-1:1
        j = order(i);
        if ~keep(j)
            continue;
        end
        k = j - 1;
        while k >= 1 && peaks(j) - peaks(k) < distance
            keep(k) = false;
            k = k - 1;
        end
        k = j + 1;
        while k <= n && peaks(k) - peaks(j) < distance
            keep(k) = false;
            k = k + 1;
        end
    end
end

function prominences = peakProminences(x, peaks)
%PEAKPROMINENCES  scipy's _peak_prominences without a window: from each peak
%   outwards while the signal stays at or below it, the lowest point on
%   either side; the prominence is the peak's height above the higher of
%   those two.
    prominences = zeros(1, numel(peaks));
    n = numel(x);
    for m = 1:numel(peaks)
        pk = peaks(m);
        higherBefore = find(x(1:pk - 1) > x(pk), 1, 'last');
        if isempty(higherBefore)
            higherBefore = 0;
        end
        higherAfter = find(x(pk + 1:n) > x(pk), 1, 'first');
        if isempty(higherAfter)
            higherAfter = n + 1;
        else
            higherAfter = pk + higherAfter;
        end
        leftMin = min(x(higherBefore + 1:pk));
        rightMin = min(x(pk:higherAfter - 1));
        prominences(m) = x(pk) - max(leftMin, rightMin);
    end
end
