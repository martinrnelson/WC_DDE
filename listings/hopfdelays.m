%HOPFDELAYS  Hopf delays of the symmetric Wilson-Cowan equilibrium
%
% For a = d, b = c, P = Q the linearisation splits into an in-phase and
% an anti-phase mode, each satisfying lambda + 1 = eta*exp(-lambda*tau).

a = -6; b = 2.5; P = 0.2; beta = 60;
f = @(z) 1./(1 + exp(-beta*z));

ustar = fzero(@(u) u - f((a+b)*u + P), [0 1]);  % equilibrium u* = v*
k  = beta*ustar*(1 - ustar);                    % slope f'(z*)
eta = k*[a+b; a-b];                              % in-phase, anti-phase
w  = sqrt(eta.^2 - 1);                           % Hopf frequencies (|eta| > 1)
j  = 0:3;
tau = (angle(eta) - atan(w) + 2*pi*j)./w;        % Hopf delays tau_j, j = 0..3

fprintf('u* = %.5f,  k = %.4f\n', ustar, k)
fprintf('in-phase:   eta = %7.3f, omega = %6.3f, tau_j =%s\n', ...
        eta(1), w(1), sprintf(' %.4f', tau(1,:)))
fprintf('anti-phase: eta = %7.3f, omega = %6.3f, tau_j =%s\n', ...
        eta(2), w(2), sprintf(' %.4f', tau(2,:)))
