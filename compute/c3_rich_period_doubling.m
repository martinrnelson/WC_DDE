%% C3  Rich parameter set: period doubling, branch switching and folds
%
% Paper:    Section 5; Figures 2(b,c), 3(b) and 4(b,c). Step 1 is Listing 5.
% Produces: data/rich_period_doubling.mat   plain arrays, for the figure scripts
%           data/rich_period_branch.mat     the DDE-BIFTOOL branch, used by c4
% Needs:    DDE-BIFTOOL (run wc_setup first)
% Run time: about 10 minutes (most of it in step 4)
%
% Steps
%   1. Find an asymmetric periodic orbit at tau = 0.28 by simulation and
%      turn it into a DDE-BIFTOOL periodic orbit.
%   2. Continue it to smaller tau and compute its Floquet multipliers.
%   3. Locate the period doubling, where a multiplier passes through -1.
%   4. Switch to the period-doubled branch and continue it.
%   5. Find the folds and further period doublings on that branch.
%   6. Save.
%
% Why simulation in step 1? This orbit does not lie on the branch that
% starts at the Hopf point (that branch carries the symmetric, anti-phase
% orbits), so we have to find it some other way.

clear; close all
p     = wc_parameters('rich');
ix    = wc_index();
funcs = wc_funcs();
mthp  = df_mthod(funcs, 'psol');     % default numerical settings for periodic orbits

%% 1. A periodic orbit from a simulation
tau0 = 0.28;
[t, u, v] = wc_simulate(p, tau0, 300, [0.3; 0.1]);

% One period is the time between successive occurrences of the largest
% maximum of u (the orbit has two maxima per period, of different heights).
late = t > 250;
tl   = t(late);
[pk, ipk] = local_maxima(u(late));
big  = ipk(pk > max(pk) - 1e-3);
T0   = mean(diff(tl(big)));
fprintf('Simulation at tau = %.2f: period %.4f\n', tau0, T0);

% A DDE-BIFTOOL periodic orbit stores one period of the solution at the
% points of a mesh on [0, 1] (time divided by the period).
deg  = 4;  ntst = 150;
mesh = linspace(0, 1, ntst*deg + 1);
times = tl(big(1)) + mesh*T0;
orbit = struct('kind', 'psol', 'parameter', wc_parvec(p, tau0), ...
    'mesh', mesh, 'degree', deg, 'period', T0, ...
    'profile', [interp1(t, u, times); interp1(t, v, times)]);

% Newton's method with tau fixed (no free parameters: the empty []) turns
% the simulated orbit into an exact periodic orbit of the collocation
% equations. p_remesh then moves the mesh points to where the orbit
% changes fastest, and we correct once more.
orbit = p_correc(funcs, orbit, [], [], mthp.point);
orbit = p_remesh(orbit, deg, ntst);
[orbit, ok] = p_correc(funcs, orbit, [], [], mthp.point);
if ~ok, error('Newton''s method did not converge for the simulated orbit.'); end

%% 2. Continue the orbit to smaller tau
% A branch needs two starting points; the second is the same orbit
% corrected at a slightly smaller delay, which sets the direction.
second = orbit;
second.parameter(ix.tau) = tau0 - 0.001;
second = p_correc(funcs, second, [], [], mthp.point);

brT = df_brnch(funcs, ix.tau, 'psol');           % empty branch, tau free
brT.parameter.min_bound = [ix.tau 0.265];
brT.parameter.max_bound = [ix.tau 0.34];
brT.parameter.max_step  = [ix.tau 0.001];
brT.method.point.adapt_mesh_after_correct = 1;
brT.point = [orbit, second];
brT = br_contn(funcs, brT, 100);
[nunst_T, ~, triv_T, brT.point] = GetStability(brT, 'funcs', funcs, ...
    'exclude_trivial', true);
fprintf('Period-T branch: trivial multiplier error at most %.1e\n', max(abs(triv_T)));

%% 3. The period doubling
% The orbit is stable at tau = 0.28 and loses stability as tau decreases.
% If the multiplier that leaves the unit circle does so through -1, this
% is a period doubling. We locate it by linear interpolation of
% g = (most negative real part of a multiplier) + 1, which is zero there.
ind_pd = find(diff(nunst_T) == 1, 1);
tauT   = arrayfun(@(q) q.parameter(ix.tau), brT.point);
g      = arrayfun(@(q) min(real(q.stability.mu)) + 1, brT.point);
s      = g(ind_pd) / (g(ind_pd) - g(ind_pd+1));
pd_point.tau = tauT(ind_pd) + s*(tauT(ind_pd+1) - tauT(ind_pd));
fprintf('Period doubling at tau = %.5f (most negative multiplier just after: %.3f)\n', ...
        pd_point.tau, g(ind_pd+1) - 1);

%% 4. Switch to the period-doubled branch
% DoublePsol builds a first guess for the period-doubled orbit: the orbit
% at point ind_pd traversed twice (period 2T), plus a small multiple
% ('radius') of the eigenfunction for the multiplier -1. It corrects this
% guess and returns a new branch with two points, ready for br_contn.
% Its mesh has twice as many intervals (300) as the period-T orbit.
[br2, ok] = DoublePsol(funcs, brT, ind_pd, 'radius', 0.01);
if ~ok, error('DoublePsol could not switch to the period-doubled branch.'); end

% Near tau = 0.268 this branch folds twice in a very narrow range of tau,
% so we take small steps. It is continued until it reaches tau = 0.265.
br2.parameter.min_bound = [ix.tau 0.265];
br2.parameter.max_bound = [ix.tau pd_point.tau + 1e-4];
br2.parameter.max_step  = [ix.tau 1e-4];
br2.method.point.adapt_mesh_after_correct = 1;
br2.method.point.print_residual_info = 0;
br2 = br_contn(funcs, br2, 300);
[nunst_2, ~, triv_2, br2.point] = GetStability(br2, 'funcs', funcs, ...
    'exclude_trivial', true);

% Accuracy check. If this number is not small (say above 1e-2), the
% multipliers cannot be trusted: increase the number of mesh intervals.
fprintf('Period-doubled branch: %d points, trivial multiplier error at most %.1e\n', ...
        numel(br2.point), max(abs(triv_2)));

%% 5. Folds and period doublings along the period-doubled branch
% A fold is where the branch turns back in tau. We keep the branch from
% the period doubling through its first two folds (beyond that it folds
% again and again in an even narrower range, which the figure does not
% need), and look for multipliers crossing -1 on that part.
doubled_all = branch_summary(br2, nunst_2, ix.tau);
[keep, ifold] = follow_branch(doubled_all.tau, 2);
doubled = take_rows(doubled_all, keep);

folds.tau    = doubled_all.tau(ifold);
folds.maxima = doubled_all.maxima(ifold, :);

g2 = cellfun(@(mu) min(real(mu)) + 1, doubled.mu);
ipd2 = find(sign(g2(1:end-1)) ~= sign(g2(2:end)));
pd2.tau = [];  pd2.maxima = [];
for k = numel(ipd2):-1:1             % backwards, so the insertions don't shift later indices
    i = ipd2(k);
    s = g2(i) / (g2(i) - g2(i+1));
    new = interpolate_row(doubled, i, s);
    pd2.tau    = [new.tau; pd2.tau];
    pd2.maxima = [new.maxima; pd2.maxima];
    % insert the bifurcation point into the branch, with the stability of
    % the side it leads into, so the solid line ends exactly there
    doubled = insert_row(doubled, i, new);
end
for k = 1:numel(folds.tau)
    fprintf('  fold of the period-doubled orbit at tau = %.5f\n', folds.tau(k));
end
for k = 1:numel(pd2.tau)
    fprintf('  further period doubling at tau = %.5f\n', pd2.tau(k));
end

%% 6. Save
period_T = branch_summary(brT, nunst_T, ix.tau);
pd_point.maxima = interp1(period_T.tau, period_T.maxima, pd_point.tau);
pd_point.mu_after  = brT.point(ind_pd+1).stability.mu;     % for Figure 3(b)
pd_point.tau_after = tauT(ind_pd+1);

% orbits for Figure 4: the asymmetric orbit at tau = 0.28, and the
% period-doubled orbit at tau = 0.272 on the stable part of its branch
profile_asym = orbit_profile(brT.point(1), ix.tau);
stable_part  = keep(1:ifold(1));
[~, k] = min(abs(doubled_all.tau(stable_part) - 0.272));
profile_doubled = orbit_profile(br2.point(stable_part(k)), ix.tau);

save(wc_file('data', 'rich_period_doubling.mat'), 'p', 'period_T', 'pd_point', ...
     'doubled', 'folds', 'pd2', 'profile_asym', 'profile_doubled');
period_T_branch = brT;               % kept for the two-parameter continuation (c4)
save(wc_file('data', 'rich_period_branch.mat'), 'period_T_branch', 'ind_pd');
fprintf('Saved data/rich_period_doubling.mat and data/rich_period_branch.mat\n');

%% ------------------------------------------------------------------------
%  Small helpers for the arrays made by branch_summary
function s = take_rows(s, idx)
% Keep only the points idx of every field.
for f = fieldnames(s)'
    s.(f{1}) = s.(f{1})(idx, :);
end
end

function new = interpolate_row(s, i, frac)
% A point a fraction frac of the way from point i to point i+1.
new.tau    = s.tau(i)    + frac*(s.tau(i+1)    - s.tau(i));
new.period = s.period(i) + frac*(s.period(i+1) - s.period(i));
new.umax   = s.umax(i)   + frac*(s.umax(i+1)   - s.umax(i));
new.umin   = s.umin(i)   + frac*(s.umin(i+1)   - s.umin(i));
new.maxima = s.maxima(i,:) + frac*(s.maxima(i+1,:) - s.maxima(i,:));
new.stable = s.stable(i+1);
new.mu     = {[]};
end

function s = insert_row(s, i, new)
% Insert the point new between points i and i+1.
for f = fieldnames(s)'
    name = f{1};
    s.(name) = [s.(name)(1:i, :); new.(name); s.(name)(i+1:end, :)];
end
end

function pr = orbit_profile(pt, itau)
% One period of an orbit, as plain arrays.
pr = struct('tau', pt.parameter(itau), 'period', pt.period, 'mesh', pt.mesh, ...
            'u', pt.profile(1,:), 'v', pt.profile(2,:));
end
