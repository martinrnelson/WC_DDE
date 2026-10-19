function [pk, idx] = local_maxima(u)
%LOCAL_MAXIMA  Local maxima of a sampled signal, e.g. from a simulation.
%
%   [pk, idx] = local_maxima(u) returns the values pk = u(idx) at the
%   points where u is larger than the sample before it and at least as
%   large as the sample after it.
%
%   For one period of a periodic orbit computed by DDE-BIFTOOL, use
%   orbit_maxima instead.

u   = u(:)';
idx = find(u(2:end-1) > u(1:end-2) & u(2:end-1) >= u(3:end)) + 1;
pk  = u(idx);
end
