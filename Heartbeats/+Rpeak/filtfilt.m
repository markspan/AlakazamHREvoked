function y = filtfilt(b, a, x)
%FILTFILT  scipy.signal.filtfilt with its defaults, the way spectHR filters:
%   forward and backward, so without phase shift.
%   Y = Rpeak.filtfilt(B, A, X). As scipy: X is extended at each end by
%   3 * max(numel(A), numel(B)) samples, reflected about its end value
%   (padtype 'odd'); the filter starts in the steady state of a step
%   (lfilter_zi) scaled to the first sample; the filtered signal is reversed,
%   filtered again the same way, reversed back, and the extensions cut off.
%   The filter itself (lfilter) is the transposed direct form II, its sums in
%   scipy's order. MATLAB's own filtfilt reflects one sample less, so its
%   ends differ.
%
%   X must be longer than the extension (scipy refuses it otherwise; spectHR
%   then leaves the signal unfiltered).
%
%   See also RPEAK.BUTTER, RPEAK.HIGHPASS.
    b = double(reshape(b, 1, []));
    a = double(reshape(a, 1, []));
    x = double(reshape(x, 1, []));
    if a(1) ~= 1
        b = b / a(1);
        a = a / a(1);
    end
    n = max(numel(a), numel(b));
    b(end + 1:n) = 0;
    a(end + 1:n) = 0;
    edge = 3 * n;
    if numel(x) <= edge
        error('Rpeak:filtfilt', ['The signal is %d samples long, and the filter needs ' ...
            'more than %d.'], numel(x), edge);
    end
    ext = [2 * x(1) - x(edge + 1:-1:2), x, 2 * x(end) - x(end - 1:-1:end - edge)];

    zi = steadyState(b, a);
    y = lfilter(b, a, ext, zi * ext(1));
    y = lfilter(b, a, fliplr(y), zi * y(end));
    y = fliplr(y);
    y = y(edge + 1:end - edge);
end

function zi = steadyState(b, a)
%STEADYSTATE  scipy.signal.lfilter_zi: the filter's state after a long run of
%   ones, from (I - A') zi = B with A the companion matrix of a.
    n = numel(a);
    companion = [-a(2:end); eye(n - 2, n - 1)];
    zi = ((eye(n - 1) - companion.') \ (b(2:end) - a(2:end) * b(1)).').';
end

function y = lfilter(b, a, x, z)
%LFILTER  scipy.signal.lfilter with the initial state Z: per sample, the
%   output from the first state, then each state from the next one, the
%   input and the output, in scipy's order of operations.
    n = numel(b);
    y = zeros(size(x));
    for k = 1:numel(x)
        xk = x(k);
        yk = z(1) + b(1) * xk;
        for m = 1:n - 2
            z(m) = z(m + 1) + xk * b(m + 1) - yk * a(m + 1);
        end
        z(n - 1) = xk * b(n) - yk * a(n);
        y(k) = yk;
    end
end
