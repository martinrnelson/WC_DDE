%% C4  Two-parameter diagram in the (b, tau) plane
%
% Paper:    Sections 3.1 and 5; Figure 5
% Produces: data/two_parameter.mat, used by figures/fig5_two_parameter.m
% Needs:    DDE-BIFTOOL, and the results of c2 and c3
% Run time: about 15 minutes (most of it in step 3)
%
% Steps
%   1. Hopf curves of the equilibrium, from the exact formula of Section 3.1.
%   2. Continue the torus bifurcation found in c2 as b and tau both vary.
%   3. Continue the period doubling found in c3 in the same way.
%   4. Find points on the torus curve where the rotation number is 1/2,
%      1/3 or 1/4.
%   5. Save.
%
% Throughout, c = b, so the two populations stay identical: wc_funcs(true)
% makes parameter 2 set both cross-couplings.

clear; close all
p     = wc_parameters('rich');
ix    = wc_index();
funcs = wc_funcs(true);
load(wc_file('data', 'rich_antiphase_branch.mat'), 'antiphase_branch', 'ind_torus');
load(wc_file('data', 'rich_period_branch.mat'),   'period_T_branch', 'ind_pd');

%% 1. Hopf curves (no continuation needed)
% For each b, the symmetric equilibrium solves u* = f((a+b)u* + P). The
% in-phase and anti-phase modes have mu = k(a+b) and mu = k(a-b), with
% k = beta u*(1 - u*), and each gives a family of Hopf delays.
f = @(z) 1 ./ (1 + exp(-p.beta*z));
hopf.b = linspace(0.05, 5.5, 500);
hopf.tau_inphase   = nan(2, numel(hopf.b));   % first 2 in-phase Hopf delays
hopf.tau_antiphase = nan(4, numel(hopf.b));   % first 4 anti-phase Hopf delays
for i = 1:numel(hopf.b)
    b = hopf.b(i);
    ustar = fzero(@(u) u - f((p.a + b)*u + p.P), [0 1]);
    k = p.beta * ustar * (1 - ustar);
    hopf.tau_inphase(:, i)   = hopf_delays_from_mu(k*(p.a + b), 2);
    hopf.tau_antiphase(:, i) = hopf_delays_from_mu(k*(p.a - b), 4);
end

% Double Hopf points: where an in-phase curve crosses an anti-phase curve
double_hopf.b = [];  double_hopf.tau = [];
for i = 1:size(hopf.tau_inphase, 1)
    for j = 1:size(hopf.tau_antiphase, 1)
        d = hopf.tau_inphase(i,:) - hopf.tau_antiphase(j,:);
        for q = find(d(1:end-1) .* d(2:end) < 0)
            s = d(q) / (d(q) - d(q+1));
            double_hopf.b(end+1)   = hopf.b(q) + s*(hopf.b(q+1) - hopf.b(q));
            double_hopf.tau(end+1) = hopf.tau_inphase(i,q) + ...
                s*(hopf.tau_inphase(i,q+1) - hopf.tau_inphase(i,q));
        end
    end
end

%% 2. The torus curve
% SetupTorusBifurcation turns the point where the anti-phase orbits lost
% stability (c2) into the first point of a branch of torus bifurcations,
% with b and tau as continuation parameters ('dir': vary b first). Each
% point on this branch is a periodic orbit with a pair of multipliers
% exactly on the unit circle, mu = exp(+-i pi omega); omega is computed
% along with the orbit.
[trfuncs, torus_branch] = SetupTorusBifurcation(funcs, antiphase_branch, ind_torus, ...
    'contpar', [ix.b ix.tau], 'dir', ix.b, 'step', 1e-3, ...
    'min_bound', [ix.b 0.02; ix.tau 0.005], ...
    'max_bound', [ix.b 5.5;  ix.tau 0.35], ...
    'max_step',  [ix.b 0.05; ix.tau 0.005], 'print_residual_info', 0);
torus_branch = br_contn(trfuncs, torus_branch, 200);   % one direction ...
torus_branch = br_rvers(torus_branch);
torus_branch = br_contn(trfuncs, torus_branch, 200);   % ... then the other

torus.b     = arrayfun(@(q) q.parameter(ix.b),   torus_branch.point);
torus.tau   = arrayfun(@(q) q.parameter(ix.tau), torus_branch.point);
torus.omega = trfuncs.get_comp(torus_branch.point, 'omega');

%% 3. The period-doubling curve
% The same idea, starting from the period doubling found in c3. This
% curve turns out to be a closed loop, so about 100 steps in each
% direction take it once round.
[pdfuncs, pd_branch] = SetupPeriodDoubling(funcs, period_T_branch, ind_pd, ...
    'contpar', [ix.b ix.tau], 'dir', ix.b, 'step', 1e-3, ...
    'min_bound', [ix.b 0.5; ix.tau 0.005], ...
    'max_bound', [ix.b 5.5; ix.tau 0.35], ...
    'max_step',  [ix.b 0.05; ix.tau 0.005], 'print_residual_info', 0);
pd_branch = br_contn(pdfuncs, pd_branch, 100);
pd_branch = br_rvers(pd_branch);
pd_branch = br_contn(pdfuncs, pd_branch, 100);

pd_curve.b   = arrayfun(@(q) q.parameter(ix.b),   pd_branch.point);
pd_curve.tau = arrayfun(@(q) q.parameter(ix.tau), pd_branch.point);

%% 4. Rotation numbers 1/2, 1/3 and 1/4 on the torus curve
% With mu = exp(+-i pi omega), a rotation number of 1/n (an angle of
% 2 pi/n) means omega = 2/n: omega = 1, 2/3 and 1/2. We locate these by
% linear interpolation between consecutive points.
% (Because the anti-phase orbit is symmetric, u(t) = v(t + T/2), whether
% these are strong resonances is decided by the half-period map; see the
% caption of Figure 5.)
resonance.b = [];  resonance.tau = [];  resonance.n = [];
for n = [2 3 4]
    d = torus.omega(:)' - 2/n;
    for q = find(d(1:end-1) .* d(2:end) < 0)
        s = d(q) / (d(q) - d(q+1));
        resonance.b(end+1)   = torus.b(q)   + s*(torus.b(q+1)   - torus.b(q));
        resonance.tau(end+1) = torus.tau(q) + s*(torus.tau(q+1) - torus.tau(q));
        resonance.n(end+1)   = n;
        fprintf('Rotation number 1/%d at b = %.4f, tau = %.5f\n', ...
                n, resonance.b(end), resonance.tau(end));
    end
end

%% 5. Save
save(wc_file('data', 'two_parameter.mat'), 'p', 'hopf', 'double_hopf', ...
     'torus', 'pd_curve', 'resonance');
fprintf('Saved data/two_parameter.mat\n');
