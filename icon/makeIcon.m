function makeIcon(out)
%MAKEICON  Draw the plugin's ribbon icon, Heartbeats.png: two beats of an ECG
%   with an event marker under each R-peak, drawn at eight times the size,
%   each stroke as the pixels within half its width of the path, and averaged
%   down to 24 x 24, in Alakazam's icon blue on a transparent ground (as the
%   Butterworth plugin draws its own).
%
%   makeIcon()      writes ../Heartbeats/Heartbeats.png
%   makeIcon(OUT)   writes OUT
    if nargin < 1
        out = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Heartbeats', 'Heartbeats.png');
    end
    scale = 8;
    n = 24 * scale;
    [x, y] = meshgrid(((1:n) - 0.5) / scale, ((1:n) - 0.5) / scale);
    % Baseline, a small P, the QRS spike, a T wave; twice.
    beat = @(x0) [x0 + 0, 13; x0 + 1.2, 13; x0 + 1.8, 12; x0 + 2.4, 13; x0 + 3.2, 13; ...
        x0 + 3.7, 14.5; x0 + 4.6, 3; x0 + 5.5, 16; x0 + 6.1, 13; x0 + 7.0, 13; ...
        x0 + 7.9, 10.8; x0 + 8.8, 13; x0 + 10, 13];
    trace = polyline([beat(1.5); beat(11.5)]);
    markers = [polyline([6.1 19; 6.1 22.5]), polyline([16.1 19; 16.1 22.5])];
    alpha = max(cat(3, stroke(x, y, trace, 1.6), stroke(x, y, markers, 1.6)), [], 3);
    alpha = squeeze(mean(mean(reshape(alpha, scale, 24, scale, 24), 1), 3));
    blue = [0x4a 0x7f 0xc9];
    rgb = repmat(reshape(uint8(blue), 1, 1, 3), 24, 24);
    imwrite(rgb, out, 'Alpha', alpha);
end

function segments = polyline(points)
%POLYLINE  Consecutive points as segments, one per column: [x1 y1 x2 y2]'.
    segments = [points(1:end - 1, :), points(2:end, :)].';
end

function a = stroke(x, y, segments, width)
%STROKE  1 within half WIDTH of any segment, 0 elsewhere (round caps).
    d = inf(size(x));
    for k = 1:size(segments, 2)
        p = segments(1:2, k).';
        q = segments(3:4, k).';
        v = q - p;
        t = ((x - p(1)) * v(1) + (y - p(2)) * v(2)) / max(sum(v .^ 2), eps);
        t = min(max(t, 0), 1);
        d = min(d, hypot(x - (p(1) + t * v(1)), y - (p(2) + t * v(2))));
    end
    a = double(d <= width / 2);
end
