function a = orbit_asymmetry(pt)
%ORBIT_ASYMMETRY  How far a periodic orbit is from being symmetric under u <-> v.
%
%   a = orbit_asymmetry(pt), where pt is a periodic orbit from DDE-BIFTOOL.
%
%   When the two populations are identical (a = d, b = c, P = Q), swapping
%   u and v maps solutions to solutions. A periodic orbit is symmetric if
%   swapping gives back the same orbit, perhaps shifted in time:
%
%       v(t) = u(t + s)   for some shift s,
%
%   with s = 0 for an in-phase orbit (u = v) and s = T/2 for an anti-phase
%   orbit. We measure
%
%       a = min over s of  max over t |v(t) - u(t + s)|,  divided by the range of u,
%
%   so a = 0 for a symmetric orbit and a > 0 for an asymmetric one. An
%   asymmetric orbit always has a mirror image (u and v exchanged) with
%   the same value of a.

n = 400;                                     % uniform time grid on one period
t = (0:n-1) / n;
[mesh, keep] = unique(pt.mesh);              % the mesh may be non-uniform
u = interp1(mesh, pt.profile(1, keep), t);
v = interp1(mesh, pt.profile(2, keep), t);

best = inf;
for k = 0:n-1
    best = min(best, max(abs(v - circshift(u, -k))));   % circshift(u,-k) is u(t + k/n)
end
a = best / (max(u) - min(u));
end
