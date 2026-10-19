%PSOL_FROM_SIM  Start a periodic-orbit computation from a simulation
%
% Requires DDE-BIFTOOL and funcs, itau from ddebiftool_hopf.m.
% Uses the beta = 60 parameters, where this orbit is not on the branch
% from the Hopf point, so it must be found another way.

par = [-6 2.5 2.5 -6 0.2 0.2 60 0.28];
fz  = @(z) 1./(1 + exp(-60*z));
dde = @(t,y,Z) [-y(1) + fz(-6*Z(1) + 2.5*Z(2) + 0.2);
                -y(2) + fz(2.5*Z(1) - 6*Z(2) + 0.2)];
sol = dde23(dde, par(itau), [0.3; 0.1], [0 300], ...
            ddeset('RelTol', 1e-6, 'AbsTol', 1e-9));

% one period = time between successive occurrences of the largest maximum
t = linspace(250, 300, 50001);  u = deval(sol, t, 1);
ip  = find(u(2:end-1) > u(1:end-2) & u(2:end-1) >= u(3:end)) + 1;
big = ip(u(ip) > max(u(ip)) - 1e-3);
T0  = t(big(2)) - t(big(1));

deg = 4; ntst = 120; mesh = linspace(0, 1, ntst*deg + 1);
po = struct('kind', 'psol', 'parameter', par, 'mesh', mesh, 'degree', deg, ...
            'profile', deval(sol, t(big(1)) + mesh*T0), 'period', T0);
mth = df_mthod(funcs, 'psol');
po  = p_correc(funcs, po, [], [], mth.point);       % Newton, tau fixed
po  = p_remesh(po, deg, ntst);                      % adapt mesh to profile
[po, suc] = p_correc(funcs, po, [], [], mth.point);

po.stability = p_stabil(funcs, po, mth.stability);
mu = po.stability.mu;  [~, i1] = min(abs(mu - 1));
fprintf('converged = %d, period = %.4f\n', suc, po.period)
fprintf('trivial multiplier = %.4f, most negative = %.4f\n', ...
        real(mu(i1)), min(real(mu)))
