function g = skew(x)
%SKEW  scipy.stats.skew(X) with its default bias=True: the third central
%   moment over the second to the power 1.5, the moments as scipy forms them
%   (the deviations squared, then times the deviations again), every mean
%   pairwise (Rpeak.npMean). NaN for a signal without spread, as scipy.
%
%   See also RPEAK.POLARITY.
    x = double(reshape(x, 1, []));
    mu = Rpeak.npMean(x);
    d = x - mu;
    d2 = d .* d;
    m2 = Rpeak.npMean(d2);
    m3 = Rpeak.npMean(d2 .* d);
    if m2 <= (eps * mu) ^ 2
        g = NaN;
    else
        g = m3 / m2 ^ 1.5;
    end
end
