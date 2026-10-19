%% C2  Rich parameter set: Hopf point, anti-phase orbits and the torus bifurcation
%
% Paper:    Section 5; Figures 2(a), 3(a) and 4(a)
% Produces: data/rich_hopf_torus.mat   plain arrays, used by the figure scripts
%           data/rich_antiphase_branch.mat   the DDE-BIFTOOL branch, used by c4
% Needs:    DDE-BIFTOOL (run wc_setup first)
% Run time: a few minutes
%
% Steps
%   1. Continue the equilibrium in tau and compute its stability.
%   2. Locate the Hopf bifurcation and compare it with Section 3.1.
%   3. Continue the anti-phase periodic orbits born there.
%   4. Find where these orbits lose stability: a torus bifurcation.
%   5. Save.
%
% Compared with c1, the sigmoid is much steeper (beta = 60), so the
% oscillations are small and fast. The steps in tau must be small, the
% first orbit must be very small ('radius'), and the mesh is adapted to
% the shape of each orbit ('adapt_mesh_after_correct').

clear; close all
p     = wc_parameters('rich');
ix    = wc_index();
funcs = wc_funcs();

%% 1. Branch of equilibria
xstar = wc_equilibrium(p);           % u* = v* = 0.0695: the populations are identical
eq = SetupStst(funcs, 'x', xstar, 'parameter', wc_parvec(p, 0.02), ...
    'contpar', ix.tau, 'step', 0.001, ...
    'min_bound', [ix.tau 0], 'max_bound', [ix.tau 0.1], ...
    'max_step', [ix.tau 0.002]);
eq = br_contn(funcs, eq, 100);
[nunst_eq, ~, ~, eq.point] = GetStability(eq, 'funcs', funcs);

%% 2. The Hopf bifurcation
ih   = find(diff(nunst_eq) == 2, 1);
hopf = p_tohopf(funcs, eq.point(ih));
mth  = df_mthod(funcs, 'hopf');
hopf = p_correc(funcs, hopf, ix.tau, [], mth.point);

H = wc_hopf_delays(p);
fprintf('Hopf point (DDE-BIFTOOL): tau = %.5f, omega = %.3f\n', ...
        hopf.parameter(ix.tau), hopf.omega);
fprintf('Hopf point (Section 3.1): tau = %.5f, omega = %.3f\n', H.tau_c, H.omega_c);

%% 3. Anti-phase periodic orbits
% With beta = 60 a change of 0.01 in u changes the argument of f by about
% 5, so the default starting radius of 0.01 is far too large: the first
% guess would not be close to an orbit and Newton's method would fail.
% We start from an orbit of size 1e-4 instead.
eq_hopf = eq;
eq_hopf.point = hopf;
po = SetupPsol(funcs, eq_hopf, 1, 'radius', 1e-4, ...
    'degree', 4, 'intervals', 40, ...
    'max_bound', [ix.tau 0.075], 'max_step', [ix.tau 5e-4]);
po.method.point.adapt_mesh_after_correct = 1;   % move mesh points to where the orbit changes fastest
po = br_contn(funcs, po, 200);

po = br_stabl(funcs, po, 0, 1);
[nunst_po, ~, triv_defect, po.point] = GetStability(po, 'funcs', funcs, ...
    'exclude_trivial', true);
fprintf('Anti-phase orbits: trivial multiplier error at most %.1e\n', ...
        max(abs(triv_defect)));

%% 4. The torus bifurcation
% The orbits are stable just after the Hopf point and later become
% unstable. ind_torus is the last stable point before that happens.
ind_torus = find(nunst_po(1:end-1) == 0 & nunst_po(2:end) > 0, 1);

% Between points ind_torus and ind_torus+1 the largest non-trivial
% multiplier crosses |mu| = 1. We locate the crossing by linear
% interpolation. A complex pair crossing means a torus bifurcation.
mu1  = nontrivial_multipliers(po.point(ind_torus).stability.mu);
mu2  = nontrivial_multipliers(po.point(ind_torus+1).stability.mu);
tau1 = po.point(ind_torus).parameter(ix.tau);
tau2 = po.point(ind_torus+1).parameter(ix.tau);
s    = (1 - max(abs(mu1))) / (max(abs(mu2)) - max(abs(mu1)));
torus_point.tau = tau1 + s*(tau2 - tau1);
torus_point.mu_after = po.point(ind_torus+1).stability.mu;   % for Figure 3(a)
torus_point.tau_after = tau2;
[~, k] = max(abs(mu2));
fprintf('Torus bifurcation at tau = %.4f; critical multipliers %s\n', ...
        torus_point.tau, num2str(mu2(k), 3));
fprintf('  angle of the critical pair: %.1f degrees\n', abs(angle(mu2(k)))*180/pi);

%% 5. Save
equilibrium = branch_summary(eq, nunst_eq, ix.tau);
orbits      = branch_summary(po, nunst_po, ix.tau);
torus_point.umax = interp1([tau1 tau2], orbits.umax(ind_torus:ind_torus+1), ...
                           torus_point.tau);
hopf_point  = struct('tau', hopf.parameter(ix.tau), 'omega', hopf.omega, ...
                     'u', hopf.x(1), 'v', hopf.x(2));

% one anti-phase orbit for Figure 4(a)
[~, k] = min(abs(orbits.tau - 0.055));
profile_antiphase = struct('tau', orbits.tau(k), 'period', po.point(k).period, ...
    'mesh', po.point(k).mesh, 'u', po.point(k).profile(1,:), ...
    'v', po.point(k).profile(2,:));

save(wc_file('data', 'rich_hopf_torus.mat'), 'p', 'equilibrium', 'orbits', ...
     'hopf_point', 'torus_point', 'profile_antiphase', 'H');
antiphase_branch = po;              % kept for the two-parameter continuation (c4)
save(wc_file('data', 'rich_antiphase_branch.mat'), 'antiphase_branch', 'ind_torus');
fprintf('Saved data/rich_hopf_torus.mat and data/rich_antiphase_branch.mat\n');

%% ------------------------------------------------------------------------
function mu = nontrivial_multipliers(mu)
% Remove the trivial multiplier, the one closest to 1.
[~, i] = min(abs(mu - 1));
mu(i) = [];
end
