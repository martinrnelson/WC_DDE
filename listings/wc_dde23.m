%WC_DDE23  Delay-induced Hopf bifurcation in the Wilson-Cowan model
%
% Simulates the model with dde23 for three delays, starting from a small
% perturbation of the equilibrium (which does not depend on tau).

a = 10; b = -10; c = 10; d = 2; P = -2; Q = -4; beta = 1;
f   = @(z) 1./(1 + exp(-beta*z));
rhs = @(t,y,Z) [-y(1) + f(a*Z(1) + b*Z(2) + P);      % Z(:,1) = y(t - tau)
                -y(2) + f(c*Z(1) + d*Z(2) + Q)];
ystar = [0.18644; 0.13377];                          % equilibrium
opts  = ddeset('RelTol', 1e-6, 'AbsTol', 1e-9);

taus = [0.05 0.12 0.5];
for k = 1:3
    sol = dde23(rhs, taus(k), ystar + [0.01; 0], [0 400], opts);
    subplot(3,1,k), plot(sol.x, sol.y(1,:)), ylabel('u')
    late = sol.y(1, sol.x > 350);
    fprintf('tau = %.2f: oscillation amplitude of u = %.4f\n', ...
            taus(k), max(late) - min(late))
end
xlabel('t')
