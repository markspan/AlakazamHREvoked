function m = npMean(x)
%NPMEAN  numpy.mean of X: its pairwise sum (Rpeak.npSum) over its length.
%
%   See also RPEAK.NPSUM.
    m = Rpeak.npSum(x) / numel(x);
end
