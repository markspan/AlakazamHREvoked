function s = npSum(x)
%NPSUM  The sum of X as numpy adds a float64 array, to the last bit.
%   S = Rpeak.npSum(X) adds the elements of X pairwise, as numpy.sum does:
%   halves until a block holds 128 or fewer, each block summed with eight
%   running sums. Bit for bit what numpy.sum gives (checked against numpy
%   2.4 on lengths from 1 to 300000), and so what numpy.mean, numpy.std and
%   numpy.median of two numbers give.
%
%   WHY IT MATTERS. spectHR estimates the sampling rate from the mean time
%   step and truncates the refractory distance in samples, int(0.3 * fs): a
%   rate an ulp below 1000 Hz gives 299 samples, an ulp above it 300. MATLAB's
%   own sum adds in another order, so it can land on the other side.
%
%   See also RPEAK.DETECT, RPEAK.NPMEDIAN, RPEAK.NPSTD.
    x = double(reshape(x, 1, []));
    s = pairwise(x, 1, numel(x));
end

function res = pairwise(x, lo, n)
%PAIRWISE  numpy's pairwise_sum on X(LO:LO+N-1).
    if n < 8
        res = 0;
        for i = 0:n - 1
            res = res + x(lo + i);
        end
    elseif n <= 128
        r = x(lo:lo + 7);
        whole = n - mod(n, 8);
        for i = 8:8:whole - 1
            r = r + x(lo + i:lo + i + 7);
        end
        res = ((r(1) + r(2)) + (r(3) + r(4))) + ((r(5) + r(6)) + (r(7) + r(8)));
        for i = whole:n - 1
            res = res + x(lo + i);
        end
    else
        half = floor(n / 2);
        half = half - mod(half, 8);
        res = pairwise(x, lo, half) + pairwise(x, lo + half, n - half);
    end
end
