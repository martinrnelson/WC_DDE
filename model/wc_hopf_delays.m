function H = wc_hopf_delays(p, nj)
%WC_HOPF_DELAYS  Delays at which the equilibrium has a Hopf bifurcation.
%
%   H = wc_hopf_delays(p)       first Hopf delay for each eigenvalue of M
%   H = wc_hopf_delays(p, nj)   the first nj Hopf delays tau_0, ..., tau_{nj-1}
%
%   Linearising the model about its equilibrium (Section 3.1) gives
%
%       delta'(t) = -delta(t) + M delta(t - tau),   M = diag(f'(z*)) W,
%
%   where W = [a b; c d] and z* is the argument of f at the equilibrium.
%   For each eigenvalue mu of M the characteristic equation is
%   lambda + 1 = mu exp(-lambda tau). A Hopf bifurcation has
%   lambda = i omega; comparing modulus and argument gives
%
%       omega = sqrt(|mu|^2 - 1),
%       omega tau = arg(mu) - arctan(omega) + 2 pi j,   j = 0, 1, 2, ...
%
%   Output fields:
%     H.xstar    the equilibrium [u*; v*]
%     H.M        the matrix M
%     H.mu       the eigenvalues of M (a column)
%     H.omega    the Hopf frequency for each eigenvalue (NaN if |mu| <= 1)
%     H.tau      row k holds the Hopf delays for eigenvalue k
%     H.tau_c    the smallest Hopf delay: where stability is first lost
%     H.omega_c  the frequency there; the orbit born there has period
%                2*pi/H.omega_c

if nargin < 2
    nj = 1;
end
W = [p.a, p.b; p.c, p.d];
I = [p.P; p.Q];
f = @(z) 1 ./ (1 + exp(-p.beta*z));

xstar = wc_equilibrium(p);
F     = f(W*xstar + I);                 % equals xstar at the equilibrium
M     = diag(p.beta * F .* (1 - F)) * W;
mu    = eig(M);

omega = nan(numel(mu), 1);
tau   = nan(numel(mu), nj);
for k = 1:numel(mu)
    [tau(k,:), omega(k)] = hopf_delays_from_mu(mu(k), nj);
end

H.xstar = xstar;
H.M     = M;
H.mu    = mu;
H.omega = omega;
H.tau   = tau;
[H.tau_c, kc] = min(tau(:,1));
H.omega_c     = omega(kc);
end
