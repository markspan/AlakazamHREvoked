function [a, b, residuals] = robustLine(x, y)
%ROBUSTLINE  The straight line y = a + b*x fitted robustly: least squares
%   with Huber's loss, by iteratively reweighted least squares, so a few
%   pairs far off the line count linearly rather than quadratically (as
%   load_xdf fits a stream's clock offsets).
%
%   [A, B, RESIDUALS] = Sync.robustLine(X, Y). The Huber threshold is 1.345
%   robust standard deviations (1.4826 times the median absolute residual),
%   the usual choice; X is centred while fitting, for precision on clock
%   times of thousands of seconds. RESIDUALS, Y less the line, are shaped
%   as Y.
%
%   See also SYNC.MATCHTRIGGERS.
    shape = size(y);
    x = double(x(:));
    y = double(y(:));
    centre = mean(x);
    design = [ones(size(x)), x - centre];
    coef = design \ y;
    for iteration = 1:50
        r = y - design * coef;
        scale = 1.4826 * median(abs(r - median(r)));
        if ~(scale > 0)
            break;
        end
        k = 1.345 * scale;
        w = sqrt(min(1, k ./ max(abs(r), eps)));   % Huber's weights, as row scales
        previous = coef;
        coef = (design .* w) \ (y .* w);
        if max(abs(coef - previous)) < 1e-12 * max(1, max(abs(coef)))
            break;
        end
    end
    b = coef(2);
    a = coef(1) - b * centre;
    residuals = reshape(y - (a + b * x), shape);
end
