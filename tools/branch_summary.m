function s = branch_summary(branch, nunst, itau)
%BRANCH_SUMMARY  Plain arrays from a DDE-BIFTOOL branch, ready for plotting.
%
%   s = branch_summary(branch, nunst, itau)
%
%   branch  a branch computed by DDE-BIFTOOL (equilibria or periodic orbits)
%   nunst   number of unstable eigenvalues/multipliers at each point, from
%           GetStability
%   itau    position of tau in the parameter vector (wc_index().tau)
%
%   For every branch:
%     s.tau     the delay at each point (a column)
%     s.stable  true where nunst == 0
%   For a branch of equilibria, also:
%     s.u, s.v  the equilibrium values
%   For a branch of periodic orbits, also:
%     s.period  the period of each orbit
%     s.umax, s.umin   the largest and smallest value of u on each orbit
%     s.maxima  the local maxima of u over one period, one column per
%               maximum, largest first. The number of columns is the most
%               common number of maxima along the branch; orbits with fewer
%               maxima get NaN.
%     s.mu      the Floquet multipliers of each orbit (a cell array)
%
%   The saved data files contain only these arrays, so the figure scripts
%   can be run without DDE-BIFTOOL.

pts      = branch.point;
s.tau    = arrayfun(@(q) q.parameter(itau), pts(:));
s.stable = nunst(:) == 0;

if strcmp(pts(1).kind, 'stst')
    s.u = arrayfun(@(q) q.x(1), pts(:));
    s.v = arrayfun(@(q) q.x(2), pts(:));
    return
end

n = numel(pts);
s.period = arrayfun(@(q) q.period, pts(:));
s.umax   = arrayfun(@(q) max(q.profile(1,:)), pts(:));
s.umin   = arrayfun(@(q) min(q.profile(1,:)), pts(:));
s.mu     = arrayfun(@(q) q.stability.mu, pts(:), 'UniformOutput', false);

peaks = cell(n, 1);  prom = cell(n, 1);
for i = 1:n
    [peaks{i}, prom{i}] = orbit_maxima(pts(i).profile(1,:));
end
count    = cellfun(@numel, peaks);
k        = mode(count);
s.maxima = nan(n, k);
for i = 1:n
    if count(i) >= k
        % keep the k most prominent peaks, then list them largest first
        [~, order] = sort(prom{i}, 'descend');
        s.maxima(i, :) = sort(peaks{i}(order(1:k)), 'descend');
    end
end
end
