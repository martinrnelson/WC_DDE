%% C3  Rich parameter set: period doubling, branch switching and folds
%
% Paper:    Section 5; Figures 2(b,c), 3(b) and 4(b,c). Step 1 is Listing 5.
% Produces: data/rich_period_doubling.mat   plain arrays, for the figure scripts
%           data/rich_period_branch.mat     the DDE-BIFTOOL branch, used by c4
% Needs:    DDE-BIFTOOL (run wc_setup first)
% Run time: about 35 minutes (most of it in steps 4, 5 and 5b)
%
% Steps
%   1. Find an asymmetric periodic orbit at tau = 0.28 by simulation and
%      turn it into a DDE-BIFTOOL periodic orbit.
%   2. Continue it to smaller tau and compute its Floquet multipliers.
%   3. Locate the period doubling, where a multiplier passes through -1.
%   4. Where does this branch come from? Follow it all the way round: it
%      is a closed loop (an isola). Watch its symmetry and stability.
%   5. Switch to the period-doubled branch and continue it; then (5b)
%      continue it further, until it reaches the second period doubling.
%   6. Find the folds and further period doublings on that branch.
%   7. Save.
%
% Why simulation in step 1? The orbits born at the Hopf points of the
% equilibrium are symmetric: in phase, u(t) = v(t), or in anti-phase,
% v(t) = u(t + T/2). This orbit is neither (Figure 4(b)), so it is not on
% any of those branches and we have to find it some other way. Step 4
% then asks where its branch begins. (Answer: nowhere. The branch is a
% closed loop, an isola, bounded by two folds near tau = 0.179 and 0.339.)
%
% The isola and the period-doubled branches are periodic orbits computed
% by continuation, so they include unstable orbits that a simulation
% would never show.

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

%% 4. Where does the period-T branch come from? Follow it all the way round
% An asymmetric orbit has a mirror image (u and v swapped), and such pairs
% typically split off a symmetric orbit at a symmetry-breaking
% bifurcation, where a Floquet multiplier passes through +1. A branch can
% also begin at a fold, where it turns back in tau. To find out which, we
% follow a copy of the branch from step 2 (brT itself is kept for c4)
% much further, 50 points at a time, and record along the way
%   * the asymmetry of each orbit (orbit_asymmetry: 0 means symmetric),
%   * every change of stability, and the multiplier that causes it,
%   * every turning point in tau (a fold).
% The answer: the branch turns back at a fold near tau = 0.179, again at a
% fold near tau = 0.339, and then returns to where it started. It is a
% closed loop of periodic orbits, an isola, so we stop once it has come
% back below its starting delay after two folds.
isola = brT;
isola.point = brT.point(1:2);                 % tau = 0.28, then 0.279: heading down
isola.parameter.min_bound = [ix.tau 0.12];
isola.parameter.max_bound = [ix.tau 0.40];
closed = false;
while ~closed && numel(isola.point) < 1000
    n_before = numel(isola.point);
    isola = br_contn(funcs, isola, 50);
    if numel(isola.point) == n_before
        warning('Continuation of the isola stopped early.');  break
    end
    tau_i = arrayfun(@(q) q.parameter(ix.tau), isola.point);
    turn  = find(diff(sign(diff(tau_i))) ~= 0) + 1;      % folds
    if numel(turn) >= 2
        back = find(tau_i(turn(2):end) <= tau_i(2), 1);   % back past the start
        if ~isempty(back)
            isola.point = isola.point(1:turn(2) + back - 1);
            closed = true;
        end
    end
end
if closed, fprintf('\nThe period-T branch closes on itself: it is an isola.\n'); end

[nunst_I, ~, triv_I, isola.point] = GetStability(isola, 'funcs', funcs, ...
    'exclude_trivial', true);
isola_summary = branch_summary(isola, nunst_I, ix.tau);
isola_summary.asymmetry = arrayfun(@orbit_asymmetry, isola.point(:));
isola_events = branch_events(isola_summary);

fprintf('Isola: %d points, tau from %.4f to %.4f, trivial multiplier error at most %.1e\n', ...
        numel(isola.point), min(isola_summary.tau), max(isola_summary.tau), max(abs(triv_I)));
fprintf('  asymmetry between %.3f and %.3f (0 would mean a symmetric orbit)\n', ...
        min(isola_summary.asymmetry), max(isola_summary.asymmetry));
fprintf('  period between %.3f and %.3f\n', min(isola_summary.period), max(isola_summary.period));
fprintf('  folds at tau = %s\n', num2str(isola_events.fold_tau', '%.5f  '));
fprintf('  period doublings at tau = %s\n', num2str(isola_events.pd_tau', '%.5f  '));
for i = find(diff(nunst_I(:)) ~= 0)'
    fprintf('  stability changes between tau = %.5f and %.5f: %s\n', ...
            isola_summary.tau(i), isola_summary.tau(i+1), ...
            crossing_type(isola.point(i).stability.mu, isola.point(i+1).stability.mu));
end
% How to read the output: if the asymmetry fell to (nearly) zero where a
% multiplier crosses +1, the asymmetric orbits would split off a symmetric
% orbit there. Instead the asymmetry stays large, and pairs of orbits are
% created and destroyed at the two folds.

%% 5. Switch to the period-doubled branch
% DoublePsol builds a first guess for the period-doubled orbit: the orbit
% at point ind_pd traversed twice (period 2T), plus a small multiple
% ('radius') of the eigenfunction for the multiplier -1. It corrects this
% guess and returns a new branch with two points, ready for br_contn.
% Its mesh has twice as many intervals (300) as the period-T orbit.
% As in c2, the default radius (0.01) is far too large for this steep
% sigmoid: Newton's method then fails to converge. We use 1e-3.
[br2, ok] = DoublePsol(funcs, brT, ind_pd, 'radius', 1e-3);
if ~ok
    error(['DoublePsol could not switch to the period-doubled branch: ', ...
           'try a smaller ''radius''.']);
end

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

%% 5b. Follow the period-doubled orbits over the whole isola
% The branch from step 5 is enough for the enlargement in Figure 2, but the
% period-doubled orbits continue beyond tau = 0.265. We carry on from its
% last two points with larger steps. There is a second period doubling on
% the isola, near tau = 0.179 (step 4), and the branch ends there.
%
% A trap: the period-T orbit traversed twice is also a periodic orbit of
% period 2T. At a period doubling the period-doubled branch meets it, and
% a continuation can carry straight on along it, round the isola. The
% test is that an orbit traversed twice repeats itself after half its
% period: cover_distance (below) is zero for such an orbit. We stop each
% period-doubled branch where it reaches it.
n_ext = 400;                                     % points on each branch
thr   = 0.01;                                    % cover_distance threshold
ext = br2;
ext.point = br2.point(end-1:end);
ext.parameter.min_bound = [ix.tau 0.17];
ext.parameter.max_bound = [ix.tau 0.35];
ext.parameter.max_step  = [ix.tau 1e-3];
ext = contn_until_cover(funcs, ext, n_ext, thr, ix.tau);
[nunst_E, ~, triv_E, ext.point] = GetStability(ext, 'funcs', funcs, 'exclude_trivial', true);
doubled_ext = branch_summary(ext, nunst_E, ix.tau);
ext_events  = branch_events(doubled_ext);
report_branch('Period-doubled branch from tau = 0.276, continued', doubled_ext, ext_events, triv_E);

% Where does this branch end? It runs down to the second period doubling
% on the isola, near tau = 0.179, and there it merges with the period-T
% orbit traversed twice. So the period-doubled orbits form a single branch
% joining the two period doublings on the isola (a 'period-doubling
% bubble'). Switching branches at tau = 0.179 with DoublePsol gives the
% same branch again, traced the other way, or so we expected; set
% explore_lower = true to switch there anyway and see which branch it finds.
explore_lower = false;
doubled_low = [];  low_events = [];
[~, j] = min(abs(isola_summary.tau(isola_events.pd_index) - 0.179));
if explore_lower
    ind_pd_low = isola_events.pd_index(j);
    for radius = [1e-3, -1e-3, 1e-2, -1e-2]
        [br3, ok] = DoublePsol(funcs, isola, ind_pd_low, 'radius', radius);
        if ~ok, continue, end
        br3.parameter.min_bound = [ix.tau 0.17];
        br3.parameter.max_bound = [ix.tau 0.35];
        br3.parameter.max_step  = [ix.tau 1e-3];
        br3.method.point.adapt_mesh_after_correct = 1;
        fprintf('\nSwitching at tau = 0.179 with radius %g\n', radius);
        [br3, left] = contn_until_cover(funcs, br3, n_ext, thr, ix.tau);
        if ~left
            fprintf('radius %g: the branch never left the orbit traversed twice\n', radius);
            continue
        end
        [nunst_L, ~, triv_L, br3.point] = GetStability(br3, 'funcs', funcs, 'exclude_trivial', true);
        doubled_low = branch_summary(br3, nunst_L, ix.tau);
        low_events  = branch_events(doubled_low);
        report_branch(sprintf('Period-doubled branch from tau = 0.179 (radius %g)', radius), ...
                      doubled_low, low_events, triv_L);
        break
    end
end
pd_low.tau    = isola_events.pd_tau(j);
pd_low.maxima = isola_events.pd_maxima(j, :);
pd_low.umax   = isola_events.pd_umax(j);

%% 6. Folds and period doublings along the period-doubled branch
% A fold is where the branch turns back in tau. We keep the branch from
% the period doubling through its first two folds (beyond that it folds
% again and again in an even narrower range, which the figure does not
% need), and look for multipliers crossing -1 on that part.
doubled_all = branch_summary(br2, nunst_2, ix.tau);
[keep, ifold] = follow_branch(doubled_all.tau, 2);
doubled = take_rows(doubled_all, keep);

folds.tau    = doubled_all.tau(ifold);
folds.maxima = doubled_all.maxima(ifold, :);
folds.umax   = doubled_all.umax(ifold);

g2 = cellfun(@(mu) min(real(mu)) + 1, doubled.mu);
ipd2 = find(sign(g2(1:end-1)) ~= sign(g2(2:end)));
pd2.tau = [];  pd2.maxima = [];  pd2.umax = [];
for k = numel(ipd2):-1:1             % backwards, so the insertions don't shift later indices
    i = ipd2(k);
    s = g2(i) / (g2(i) - g2(i+1));
    new = interpolate_row(doubled, i, s);
    pd2.tau    = [new.tau; pd2.tau];
    pd2.maxima = [new.maxima; pd2.maxima];
    pd2.umax   = [new.umax; pd2.umax];
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

%% 7. Save
period_T = branch_summary(brT, nunst_T, ix.tau);
pd_point.maxima = interp1(period_T.tau, period_T.maxima, pd_point.tau);
pd_point.umax   = interp1(period_T.tau, period_T.umax, pd_point.tau);
pd_point.mu_after  = brT.point(ind_pd+1).stability.mu;     % for Figure 3(b)
pd_point.tau_after = tauT(ind_pd+1);

% orbits for Figure 4: the asymmetric orbit at tau = 0.28, and the
% period-doubled orbit at tau = 0.272 on the stable part of its branch
profile_asym = orbit_profile(brT.point(1), ix.tau);
stable_part  = keep(1:ifold(1));
[~, k] = min(abs(doubled_all.tau(stable_part) - 0.272));
profile_doubled = orbit_profile(br2.point(stable_part(k)), ix.tau);

save(wc_file('data', 'rich_period_doubling.mat'), 'p', 'period_T', 'pd_point', ...
     'doubled', 'folds', 'pd2', 'profile_asym', 'profile_doubled', ...
     'isola_summary', 'isola_events', 'pd_low', ...
     'doubled_ext', 'ext_events', 'doubled_low', 'low_events');
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

function txt = crossing_type(~, mu_after)
% Which multiplier has just crossed the unit circle: the non-trivial
% multiplier closest to the unit circle just after the change of stability.
mu_after = remove_trivial(mu_after);
[~, k] = min(abs(abs(mu_after) - 1));
mu = mu_after(k);
if abs(imag(mu)) > 1e-6
    txt = sprintf('complex pair %s (torus)', num2str(mu, 3));
elseif real(mu) > 0
    txt = sprintf('real multiplier %.3f near +1 (fold or symmetry-breaking)', real(mu));
else
    txt = sprintf('real multiplier %.3f near -1 (period doubling)', real(mu));
end
end

function mu = remove_trivial(mu)
[~, i] = min(abs(mu - 1));
mu(i) = [];
end

function pr = orbit_profile(pt, itau)
% One period of an orbit, as plain arrays.
pr = struct('tau', pt.parameter(itau), 'period', pt.period, 'mesh', pt.mesh, ...
            'u', pt.profile(1,:), 'v', pt.profile(2,:));
end

function ev = branch_events(s)
% Folds and period doublings along a branch summary (from branch_summary).
% A fold is a turning point in tau; a period doubling is where the most
% negative real multiplier passes through -1 (found by interpolation).
ev.fold_index = find(diff(sign(diff(s.tau))) ~= 0) + 1;
ev.fold_tau    = s.tau(ev.fold_index);
ev.fold_maxima = s.maxima(ev.fold_index, :);
ev.fold_umax   = s.umax(ev.fold_index);
g = cellfun(@(mu) min(real(mu)) + 1, s.mu);
ev.pd_index = find(sign(g(1:end-1)) ~= sign(g(2:end)));
frac = g(ev.pd_index) ./ (g(ev.pd_index) - g(ev.pd_index + 1));
ev.pd_tau    = s.tau(ev.pd_index) + frac .* (s.tau(ev.pd_index+1) - s.tau(ev.pd_index));
ev.pd_maxima = s.maxima(ev.pd_index, :) + ...
    frac .* (s.maxima(ev.pd_index+1, :) - s.maxima(ev.pd_index, :));
ev.pd_umax   = s.umax(ev.pd_index) + frac .* (s.umax(ev.pd_index+1) - s.umax(ev.pd_index));
end

function report_branch(name, s, ev, triv)
fprintf('\n%s: %d points, tau from %.4f to %.4f (range %.4f to %.4f)\n', name, ...
        numel(s.tau), s.tau(1), s.tau(end), min(s.tau), max(s.tau));
fprintf('  trivial multiplier error at most %.1e\n', max(abs(triv)));
fprintf('  %d folds', numel(ev.fold_tau));
if ~isempty(ev.fold_tau)
    fprintf(', between tau = %.5f and %.5f', min(ev.fold_tau), max(ev.fold_tau));
end
fprintf('\n  period doublings at tau = %s\n', num2str(ev.pd_tau', '%.5f  '));
fprintf('  stable at %d of its points\n', sum(s.stable));
end

function d = cover_distance(pt)
% Zero if the orbit pt repeats itself after half its period, as the
% period-T orbit traversed twice does: max |x(t) - x(t + T/2)| over the
% orbit, divided by the range of u.
n = 400;  t = (0:n-1)/n;
[mesh, keep] = unique(pt.mesh);
x = interp1(mesh, pt.profile(:, keep)', t)';
d = max(max(abs(x - circshift(x, -n/2, 2)))) / (max(x(1,:)) - min(x(1,:)));
end

function [br, left] = contn_until_cover(funcs, br, n, thr, itau)
% Continue a period-doubled branch 10 points at a time, and stop as soon
% as it reaches the period-T orbit traversed twice (cover_distance below
% thr), so that it does not carry on round the isola. A branch that starts
% at a period doubling starts on that orbit, so we first wait for it to
% move away; left is false if it never does within the first 30 points.
left = cover_distance(br.point(end)) > thr;
while numel(br.point) < n
    n_before = numel(br.point);
    br = br_contn(funcs, br, 10);
    if numel(br.point) == n_before, break, end          % continuation stopped
    d = arrayfun(@cover_distance, br.point(n_before+1:end));
    fprintf('  tau = %.5f, distance from the orbit traversed twice %.3g\n', ...
            br.point(end).parameter(itau), d(end));
    if ~left
        k = find(d > thr, 1);
        if isempty(k)
            if numel(br.point) >= 30, return, end       % never left
            continue
        end
        left = true;  d(1:k) = inf;                     % only look beyond here
    end
    back = find(d < thr, 1);
    if ~isempty(back)
        br.point = br.point(1:n_before + back);
        fprintf('  reached the period-T orbit traversed twice at tau = %.5f: stopped\n', ...
                br.point(end).parameter(itau));
        return
    end
end
end