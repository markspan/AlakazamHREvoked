function s = npStd(x)
%NPSTD  numpy.std of X: the population standard deviation (ddof 0), computed
%   as numpy computes it, the mean of the squared deviations from the mean,
%   both means pairwise (Rpeak.npMean). MATLAB's std divides by n - 1.
%
%   See also RPEAK.NPMEAN.
    x = double(reshape(x, 1, []));
    d = x - Rpeak.npMean(x);
    s = sqrt(Rpeak.npMean(d .* d));
end
