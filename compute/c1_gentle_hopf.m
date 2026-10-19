%% C1  Gentle parameter set: equilibrium, Hopf point and periodic orbits
%
% Paper:    Section 5, Figure 1. This is Listing 4 with comments added and
%           the results saved.
% Produces: data/gentle_hopf.mat, used by figures/fig1_gentle_hopf.m
% Needs:    DDE-BIFTOOL (run wc_setup first)
% Run time: about a minute
%
% Steps
%   1. Continue the equilibrium as tau increases and compute its stability.
%   2. Locate the Hopf bifurcation and compare it with Section 3.1.
%   3. Continue the branch of periodic orbits born at the Hopf point.
%   4. Save the results as plain arrays.

clear; close all
p     = wc_parameters('gentle');
ix    = wc_index();          % ix.tau = 8: position of tau in the parameter vector
funcs = wc_funcs();          % the model in DDE-BIFTOOL form

%% 1. Branch of equilibria
% SetupStst creates a branch of equilibria starting at the given point,
% with tau ('contpar') as the continuation parameter; br_contn then takes
% up to 200 steps along it. The equilibrium itself does not change with
% tau, so this branch is a horizontal line: what changes is its stability.
xstar = wc_equilibrium(p);
eq = SetupStst(funcs, 'x', xstar, 'parameter', wc_parvec(p, 0.01), ...
    'contpar', ix.tau, 'step', 0.005, ...
    'min_bound', [ix.tau 0], 'max_bound', [ix.tau 1], ...
    'max_step', [ix.tau 0.01]);
eq = br_contn(funcs, eq, 200);

% GetStability computes the rightmost roots of the characteristic
% equation at each point; nunst_eq(k) is the number with positive real part.
[nunst_eq, ~, ~, eq.point] = GetStability(eq, 'funcs', funcs);

%% 2. The Hopf bifurcation
% Stability is lost where a pair of complex roots crosses the imaginary
% axis, so nunst jumps by 2. p_tohopf turns the point just before the jump
% into a first guess for the Hopf point, and p_correc solves for it exactly.
ih   = find(diff(nunst_eq) == 2, 1);
hopf = p_tohopf(funcs, eq.point(ih));
mth  = df_mthod(funcs, 'hopf');
hopf = p_correc(funcs, hopf, ix.tau, [], mth.point);

H = wc_hopf_delays(p);       % the exact values from Section 3.1
fprintf('Hopf point (DDE-BIFTOOL): tau = %.5f, omega = %.4f\n', ...
        hopf.parameter(ix.tau), hopf.omega);
fprintf('Hopf point (Section 3.1): tau = %.5f, omega = %.4f\n', ...
        H.tau_c, H.omega_c);

%% 3. Periodic orbits born at the Hopf point
% SetupPsol starts a branch of small periodic orbits at the Hopf point.
% Each orbit is stored as piecewise polynomials of the given 'degree' on a
% mesh of 'intervals' pieces covering one period; 'radius' is the size of
% the first, tiny orbit.
eq_hopf = eq;
eq_hopf.point = hopf;        % a branch whose only point is the Hopf point
po = SetupPsol(funcs, eq_hopf, 1, 'radius', 0.01, ...
    'degree', 4, 'intervals', 30, ...
    'max_bound', [ix.tau 1], 'max_step', [ix.tau 0.01]);
po = br_contn(funcs, po, 150);

% The stability of a periodic orbit is decided by its Floquet multipliers
% mu: it is stable if all of them satisfy |mu| < 1. One multiplier is
% always exactly 1 (the 'trivial' multiplier); 'exclude_trivial' leaves it
% out of the count, and triv_defect is how far its computed value is from
% 1, which is a useful check on the accuracy of the computation.
po = br_stabl(funcs, po, 0, 1);
[nunst_po, ~, triv_defect, po.point] = GetStability(po, 'funcs', funcs, ...
    'exclude_trivial', true);
fprintf('Periodic orbits: tau from %.3f to %.3f; trivial multiplier error at most %.1e\n', ...
        po.point(1).parameter(ix.tau), po.point(end).parameter(ix.tau), ...
        max(abs(triv_defect)));

%% 4. Save plain arrays for plotting
equilibrium = branch_summary(eq, nunst_eq, ix.tau);
orbits      = branch_summary(po, nunst_po, ix.tau);
hopf_point  = struct('tau', hopf.parameter(ix.tau), 'omega', hopf.omega, ...
                     'u', hopf.x(1), 'v', hopf.x(2));
save(wc_file('data', 'gentle_hopf.mat'), 'p', 'equilibrium', 'orbits', ...
     'hopf_point', 'H');
fprintf('Saved data/gentle_hopf.mat\n');
