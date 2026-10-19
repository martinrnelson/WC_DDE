function [pk, prominence] = orbit_maxima(u)
%ORBIT_MAXIMA  Local maxima of u over one period of a periodic orbit.
%
%   [pk, prominence] = orbit_maxima(u), where u is the first row of a
%   DDE-BIFTOOL profile, p.profile(1,:).
%
%   Two differences from local_maxima:
%   * The signal is periodic, so it wraps round: the last sample is a
%     repeat of the first and is dropped, and a peak may sit at the join.
%   * Near folds the computed orbits can have tiny numerical wiggles. Each
%     peak's prominence is how far it rises above the higher of the two
%     lowest points on either side before a taller peak is reached; peaks
%     less prominent than 2% of the oscillation range are ignored.
%
%   pk is sorted from largest to smallest.

u = u(:)';
u = u(1:end-1);                          % last point = first point
n = numel(u);
up = circshift(u, 1);  un = circshift(u, -1);
candidates = find(u > up & u >= un);
range = max(u) - min(u);

pk = [];  prominence = [];
for j = candidates
    lowest = [u(j), u(j)];               % lowest point to the left, right
    for side = 1:2
        step = 2*side - 3;               % -1 (left) or +1 (right)
        k = j;
        for count = 1:n-1
            k = mod(k - 1 + step, n) + 1;
            if u(k) > u(j), break, end   % reached a taller peak
            lowest(side) = min(lowest(side), u(k));
        end
    end
    prom = u(j) - max(lowest);
    if prom > 0.02 * range
        pk(end+1) = u(j);                %#ok<AGROW>
        prominence(end+1) = prom;        %#ok<AGROW>
    end
end
[pk, order] = sort(pk, 'descend');
prominence = prominence(order);
end
