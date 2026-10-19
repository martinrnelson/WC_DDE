%DDEBIFTOOL_HOPF  Equilibria, Hopf point and periodic orbits as tau varies
%
% Requires DDE-BIFTOOL (v3.x) on the MATLAB path.
% Parameters p = [a b c d P Q beta tau]; the delay is parameter 8.

par = [10 -10 10 2 -2 -4 1 0.01];  itau = 8;
f   = @(z,beta) 1./(1 + exp(-beta*z));
rhs = @(x,p) [-x(1,1) + f(p(1)*x(1,2) + p(2)*x(2,2) + p(5), p(7));
              -x(2,1) + f(p(3)*x(1,2) + p(4)*x(2,2) + p(6), p(7))];
funcs = set_funcs('sys_rhs', rhs, 'sys_tau', @() itau);  % x(:,2) = x(t-tau)

% 1. Branch of equilibria in tau, with stability
stst = SetupStst(funcs, 'x', [0.18644; 0.13377], 'parameter', par, ...
    'contpar', itau, 'step', 0.002, 'max_bound', [itau 0.5], ...
    'max_step', [itau 0.005]);
stst = br_contn(funcs, stst, 100);
[nunst, ~, ~, stst.point] = GetStability(stst, 'funcs', funcs);

% 2. Locate the Hopf point, where two roots cross the imaginary axis
ih   = find(diff(nunst) == 2, 1);
hopf = p_tohopf(funcs, stst.point(ih));
mth  = df_mthod(funcs, 'hopf');
hopf = p_correc(funcs, hopf, itau, [], mth.point);
fprintf('Hopf at tau = %.5f, omega = %.4f\n', hopf.parameter(itau), hopf.omega)

% 3. Branch of periodic orbits from the Hopf point, with Floquet multipliers
po = SetupPsol(funcs, stst, ih, 'degree', 4, 'intervals', 30, ...
    'max_bound', [itau 1], 'max_step', [itau 0.01]);
po = br_contn(funcs, po, 150);
po = br_stabl(funcs, po, 0, 1);
[nunst_po, ~, ~, po.point] = GetStability(po, 'funcs', funcs, ...
                                          'exclude_trivial', true);
tau = arrayfun(@(p) p.parameter(itau), po.point);
amp = arrayfun(@(p) max(p.profile(1,:)) - min(p.profile(1,:)), po.point);
plot(tau, amp, 'b.-'), xlabel('\tau'), ylabel('amplitude of u')
fprintf('unstable periodic orbits found: %d\n', nnz(nunst_po > 0))
