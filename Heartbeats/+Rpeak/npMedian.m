function m = npMedian(x)
%NPMEDIAN  numpy.median of X: the middle value, or the mean of the two middle
%   values, (a + b) / 2, as numpy takes it. NaN when X holds a NaN.
%
%   See also RPEAK.NPMEAN.
    x = sort(double(reshape(x, 1, [])));
    n = numel(x);
    if n == 0 || isnan(x(end))
        m = NaN;
    elseif mod(n, 2)
        m = x((n + 1) / 2);
    else
        m = (x(n / 2) + x(n / 2 + 1)) / 2;
    end
end
