function m = matchTriggers(tHere, codesHere, tThere, codesThere, shared, tolerance)
%MATCHTRIGGERS  Pair the triggers two recordings both received, and the line
%   that carries one recording's clock onto the other's.
%
%   M = Sync.matchTriggers(HERE, CODESHERE, THERE, CODESTHERE, SHARED, TOLERANCE)
%   where "here" is the recording times are carried to (the EEG) and "there"
%   the one they come from (the ECG): HERE and THERE are trigger times in
%   seconds, each on its recording's own clock, CODES* their codes as
%   Sync.normaliseCode writes them, SHARED the codes to use ({} for every
%   code both have) and TOLERANCE the furthest a pair may lie from the line
%   at first (seconds).
%
%   The line is tHere = a + b * tThere: a the offset between the clocks, b
%   their rate ratio (two free-running quartz clocks differ by tens of parts
%   per million, tens of ms an hour). Found in four steps, so that triggers
%   missing from either recording, or extra in one, do no harm:
%     1. A VOTE. Every pair of same-code triggers, one from each recording,
%        proposes an offset for each rate ratio tried; the true line gathers
%        nearly every right pair in one narrow bin, wrong pairs scatter. A
%        coarse grid of ratios (1 +- 300 ppm, 20 ppm steps, 50 ms bins), then
%        a fine one around the best (2 ppm steps, 10 ms bins).
%     2. PAIRS. Each trigger here is paired with the same-code trigger there
%        that the line puts nearest it, within the tolerance; a trigger there
%        claimed twice goes to the nearer.
%     3. THE LINE fitted to the pairs, robustly (Sync.robustLine).
%     4. Again from 2 with the new line, the tolerance narrowed to six
%        robust standard deviations of the misfit (2 ms at least), until the
%        pairs no longer change.
%
%   M has a, b, the pairs (indices pairsHere and pairsThere), their
%   residuals (seconds, on here's clock), shared (the codes used), ambiguity
%   (the second best vote, at another offset, over the best: near 1 when the
%   triggers are too regular to tell one offset from another) and problem
%   (why no line could be found, '' when one was).
%
%   See also SYNC.ROBUSTLINE, SYNC.FROMRECORDING.
    tHere = double(reshape(tHere, 1, []));
    tThere = double(reshape(tThere, 1, []));
    codesHere = reshape(cellstr(codesHere), 1, []);
    codesThere = reshape(cellstr(codesThere), 1, []);
    m = struct('a', NaN, 'b', NaN, 'pairsHere', zeros(1, 0), 'pairsThere', zeros(1, 0), ...
        'residuals', zeros(1, 0), 'shared', {{}}, 'ambiguity', NaN, 'problem', '');

    common = intersect(unique(codesHere), unique(codesThere));
    if ~isempty(shared)
        common = intersect(common, reshape(cellstr(shared), 1, []));
    end
    m.shared = reshape(common, 1, []);
    if isempty(common)
        m.problem = 'the two recordings have no trigger code in common';
        return;
    end

    % Every same-code pair, for the vote.
    allHere = zeros(0, 1);
    allThere = zeros(0, 1);
    for c = m.shared
        here = find(strcmp(codesHere, c{1}));
        there = find(strcmp(codesThere, c{1}));
        [I, J] = ndgrid(here, there);
        allHere = [allHere; I(:)]; %#ok<AGROW>
        allThere = [allThere; J(:)]; %#ok<AGROW>
        if numel(allHere) > 2e7
            m.problem = 'there are too many triggers of one code to pair them all';
            return;
        end
    end
    x = reshape(tThere(allThere), [], 1);
    y = reshape(tHere(allHere), [], 1);

    ppm = 1e-6;
    ratio = vote(x, y, 1 + (-300:20:300) * ppm, 0.05);
    [ratio, offset, best, second] = vote(x, y, ratio + (-20:2:20) * ppm, 0.01);
    m.ambiguity = second / best;
    a = offset;
    b = ratio;

    pairsHere = zeros(1, 0);
    pairsThere = zeros(1, 0);
    residuals = zeros(1, 0);
    for iteration = 1:20
        [newHere, newThere] = pairUp(tHere, codesHere, tThere, codesThere, m.shared, a, b, tolerance);
        if numel(newHere) < 2
            m.problem = sprintf('only %d trigger(s) could be paired', numel(newHere));
            m.pairsHere = newHere;
            m.pairsThere = newThere;
            return;
        end
        [a, b, residuals] = Sync.robustLine(tThere(newThere), tHere(newHere));
        scale = 1.4826 * median(abs(residuals - median(residuals)));
        narrower = min(tolerance, max(0.002, 6 * scale));
        settled = isequal(newHere, pairsHere) && isequal(newThere, pairsThere) && narrower == tolerance;
        pairsHere = newHere;
        pairsThere = newThere;
        tolerance = narrower;
        if settled
            break;
        end
    end
    m.a = a;
    m.b = b;
    m.pairsHere = pairsHere;
    m.pairsThere = pairsThere;
    m.residuals = reshape(residuals, 1, []);
end

% ======================================================================= %
function [ratio, offset, best, second] = vote(x, y, ratios, bin)
%VOTE  For each rate ratio, the offset most pairs agree on (counted over two
%   adjacent bins of width BIN); the ratio whose best offset gathers most.
%   SECOND is the most votes, at that ratio, for an offset further than four
%   bins from the best.
    best = -1;
    second = 0;
    ratio = ratios(1);
    offset = 0;
    for r = ratios
        o = y - r * x;
        keys = sort(round(o / bin));
        [values, first] = unique(keys);
        counts = diff([first; numel(keys) + 1]);
        pairCounts = counts;
        next = [diff(values) == 1; false];
        pairCounts(next) = counts(next) + counts(find(next) + 1);
        [top, at] = max(pairCounts);
        if top > best
            best = top;
            ratio = r;
            centre = values(at) * bin;
            offset = median(o(o >= centre - 0.5 * bin & o < centre + 1.5 * bin));
            far = abs(values - values(at)) > 4;
            second = max([0; pairCounts(far)]);
        end
    end
end

function [pairsHere, pairsThere] = pairUp(tHere, codesHere, tThere, codesThere, common, a, b, tolerance)
%PAIRUP  Each trigger here with the same-code trigger there that the line
%   puts nearest it, within TOLERANCE; one there claimed twice keeps the
%   nearer.
    pairsHere = zeros(1, 0);
    pairsThere = zeros(1, 0);
    distance = zeros(1, 0);
    for c = common
        here = find(strcmp(codesHere, c{1}));
        there = find(strcmp(codesThere, c{1}));
        if isempty(here) || isempty(there)
            continue;
        end
        mapped = a + b * tThere(there);                 % there's triggers on here's clock
        [d, k] = min(abs(tHere(here).' - mapped), [], 2);
        ok = d.' <= tolerance;
        pairsHere = [pairsHere, here(ok)]; %#ok<AGROW>
        pairsThere = [pairsThere, there(k(ok))]; %#ok<AGROW>
        distance = [distance, d(ok).']; %#ok<AGROW>
    end
    [~, order] = sort(distance);
    [~, firstClaim] = unique(pairsThere(order), 'stable');
    keep = order(firstClaim);
    [pairsHere, byTime] = sort(pairsHere(keep));
    pairsThere = pairsThere(keep(byTime));
end
