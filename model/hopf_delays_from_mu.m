function [tau, omega] = hopf_delays_from_mu(mu, nj)
%HOPF_DELAYS_FROM_MU  Hopf delays of  lambda + 1 = mu exp(-lambda tau).
%
%   [tau, omega] = hopf_delays_from_mu(mu, nj) returns the first nj delays
%   tau(1) < tau(2) < ... at which lambda = i omega is a root, and the
%   frequency omega. If |mu| <= 1 there are none and NaN is returned.
%
%   From Section 3.1:  omega = sqrt(|mu|^2 - 1)  (modulus), and
%   omega tau = arg(mu) - arctan(omega) + 2 pi j  (argument). The mod
%   picks the smallest positive value for j = 0; the delays are then
%   equally spaced, 2 pi / omega apart.

if nargin < 2
    nj = 1;
end
if abs(mu) <= 1
    tau = nan(1, nj);  omega = NaN;
    return
end
omega = sqrt(abs(mu)^2 - 1);
j     = 0:nj-1;
tau   = (mod(angle(mu) - atan(omega), 2*pi) + 2*pi*j) / omega;
end
