%DDEBIFTOOL_HOPF_CHAOTIC  Hopf point and torus bifurcation as tau varies
%
% Requires DDE-BIFTOOL (v3.x) on the MATLAB path.
% Parameters p = [a b c d P Q beta tau] of Coombes & Laing (2009), Figs 8-9;
% the delay is parameter 8.

par = [-6 2.5 2.5 -6 0.2 0.2 60 0.02];  itau = 8;
f   = @(z,beta) 1./(1 + exp(-beta*z));
rhs = @(x,p) [-x(1,1) + f(p(1)*x(1,2) + p(2)*x(2,2) + p(5), p(7));
              -x(2,1) + f(p(3)*x(1,2) + p(4)*x(2,2) + p(6), p(7))];
funcs = set_funcs('sys_rhs', rhs, 'sys_tau', @() itau);  % x(:,2) = x(t-tau)

% 1. Branch of equilibria in tau, with stability
stst = SetupStst(funcs, 'x', [0.069497; 0.069497], 'parameter', par, ...
    'contpar', itau, 'step', 0.001, 'max_bound', [itau 0.1], ...
    'max_step', [itau 0.002]);
stst = br_contn(funcs, stst, 100);
[nunst, ~, ~, stst.point] = GetStability(stst, 'funcs', funcs);

% 2. Locate the Hopf point, where two roots cross the imaginary axis
ih   = find(diff(nunst) == 2, 1);
hopf = p_tohopf(funcs, stst.point(ih));
mth  = df_mthod(funcs, 'hopf');
hopf = p_correc(funcs, hopf, itau, [], mth.point);
fprintf('Hopf at tau = %.5f, omega = %.3f\n', hopf.parameter(itau), hopf.omega)

% 3. Periodic orbits from the Hopf point. The sigmoid is steep (beta = 60),
%    so start with a very small orbit and take small steps in tau.
stst.point = hopf;
po = SetupPsol(funcs, stst, 1, 'radius', 1e-4, 'degree', 4, ...
    'intervals', 40, 'max_bound', [itau 0.075], 'max_step', [itau 5e-4]);
po = br_contn(funcs, po, 200);
po = br_stabl(funcs, po, 0, 1);
[nunst_po, ~, ~, po.point] = GetStability(po, 'funcs', funcs, ...
                                          'exclude_trivial', true);

% 4. First loss of stability of the orbits, and the multipliers just after it
tau = arrayfun(@(p) p.parameter(itau), po.point);
it  = find(nunst_po > 0, 1);
mu  = po.point(it).stability.mu;   [~, k] = max(abs(mu));
fprintf('orbits lose stability near tau = %.4f; leading multiplier %s\n', ...
        tau(it), num2str(mu(k), 3))
amp = arrayfun(@(p) max(p.profile(1,:)) - min(p.profile(1,:)), po.point);
plot(tau, amp, 'b.-'), xlabel('\tau'), ylabel('amplitude of u')
