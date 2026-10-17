function [b, a] = butter(order, wn, kind)
%BUTTER  scipy.signal.butter for a digital filter, the way spectHR designs its
%   ECG filters.
%   [B, A] = Rpeak.butter(ORDER, WN, KIND) is a Butterworth filter of ORDER,
%   KIND 'high' (WN a scalar) or 'band' (WN = [low high]), WN a fraction of
%   the Nyquist frequency. scipy's steps, in its order: the analogue
%   prototype's poles (buttap), the frequencies prewarped, the prototype
%   moved to a high-pass (lp2hp_zpk) or band-pass (lp2bp_zpk), the bilinear
%   transform (bilinear_zpk), and the polynomials (zpk2tf). Complex numbers
%   are divided as numpy divides them (Smith's method). Needs no Signal
%   Processing Toolbox.
%
%   TO THE LAST BIT OR TWO. numpy multiplies complex numbers with fused
%   multiply-adds where the processor has them (in np.prod, and in the BLAS
%   its polynomial convolutions call), which MATLAB cannot reproduce, so a
%   coefficient can differ from scipy's in its last one or two bits, as
%   spectHR's own do between processors. On a filter as sharp as spectHR's
%   0.5 Hz high-pass that is a few parts in ten million of the signal, and
%   the R-peaks do not move (SpectHRParityTest).
%
%   See also RPEAK.FILTFILT, RPEAK.HIGHPASS, RPEAK.POLARITY.
    m = (-order + 1):2:(order - 1);
    p = -exp(1i * pi * m / (2 * order));      % buttap: no zeros, gain 1
    z = complex(zeros(1, 0));
    k = 1;

    fs = 2;
    warped = 2 * fs * tan(pi * wn / fs);
    switch kind
        case 'high'
            degree = numel(p) - numel(z);
            k = k * real(divide(prod(-z), prod(-p)));
            p = divide(complex(warped), p);
            z = complex(zeros(1, degree));
        case 'band'
            bw = warped(2) - warped(1);
            wo = sqrt(warped(1) * warped(2));
            degree = numel(p) - numel(z);
            pl = p * bw / 2;
            p = [pl + sqrt(pl .* pl - wo * wo), pl - sqrt(pl .* pl - wo * wo)];
            z = complex(zeros(1, degree));
            k = k * bw ^ degree;
        otherwise
            error('Rpeak:butter', 'KIND is ''high'' or ''band''.');
    end

    % The bilinear transform.
    fs2 = 2 * fs;
    degree = numel(p) - numel(z);
    k = k * real(divide(prod(fs2 - z), prod(fs2 - p)));
    z = [divide(fs2 + z, fs2 - z), -ones(1, degree)];
    p = divide(fs2 + p, fs2 - p);

    b = real(k * fromRoots(z));
    a = real(fromRoots(p));
end

function q = divide(x, y)
%DIVIDE  X ./ Y for complex numbers as numpy's complex division does it:
%   Smith's method, scaling by the larger part of the divisor.
    x = complex(x) + zeros(size(y));
    y = complex(y) + zeros(size(x));
    q = complex(zeros(size(y)));
    for j = 1:numel(y)
        ar = real(x(j)); ai = imag(x(j));
        br = real(y(j)); bi = imag(y(j));
        if abs(br) >= abs(bi)
            rat = bi / br;
            scl = 1 / (br + bi * rat);
            q(j) = complex((ar + ai * rat) * scl, (ai - ar * rat) * scl);
        else
            rat = br / bi;
            scl = 1 / (bi + br * rat);
            q(j) = complex((ar * rat + ai) * scl, (ai * rat - ar) * scl);
        end
    end
end

function c = fromRoots(r)
%FROMROOTS  numpy.poly: the monic polynomial with roots R, one root at a time,
%   each a convolution with [1, -root].
    c = complex(1);
    for j = 1:numel(r)
        c = [c, 0] + [0, c * -r(j)];
    end
end
