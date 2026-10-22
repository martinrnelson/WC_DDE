%% FIGURE 2  Secondary bifurcations of periodic orbits, rich parameter set
%
% (a) Hopf bifurcation of the equilibrium and torus bifurcation of the
%     anti-phase orbits; past the torus bifurcation, local maxima from the
%     simulations of c5 can be added if desired
%     (if data/route_to_chaos.mat exists)
% (b) the isola of asymmetric periodic orbits (blue), bounded by two folds,
%     with its two period doublings and the period-doubled orbits (orange)
% (b)-(d) show the maximum of u over each orbit; (a) shows local maxima
% (c) enlargement of the shaded region in (b): the period doubling at
%     tau = 0.276 and the folds of the period-doubled branch
% (d) enlargement of the shaded region in (c)
% Data: data/rich_hopf_torus.mat (c2), data/rich_period_doubling.mat (c3)
%       and, for panel (a), data/route_to_chaos.mat (c5)

clear; close all
A = load(wc_file('data', 'rich_hopf_torus.mat'));
B = load(wc_file('data', 'rich_period_doubling.mat'));
S = paper_style();
% c3 can also switch branches at the lower period doubling (explore_lower);
% that gives the same bubble traced the other way, so it is not drawn
% unless you ask for it here.
show_lower_switch = false;
if ~show_lower_switch, B.doubled_low = [];  B.low_events = []; end

fig = new_figure(1.25*S.width, 12);
% panel positions [left bottom width height], as fractions of the figure,
% with enough space between the rows for the tau labels and panel titles
pos = {[0.10 0.60 0.37 0.33], [0.59 0.60 0.37 0.33], ...
       [0.10 0.11 0.37 0.33], [0.59 0.11 0.37 0.33]};

%% (a) equilibrium, anti-phase orbits, Hopf and torus points
ax = axes(fig, 'Position', pos{1});  hold(ax, 'on');
plot_branch(ax, A.equilibrium.tau, A.equilibrium.u, A.equilibrium.stable, S.eq);
for k = 1:size(A.orbits.maxima, 2)        % one maximum per period here
    plot_branch(ax, A.orbits.tau, A.orbits.maxima(:, k), A.orbits.stable, S.po);
end
% Past the torus bifurcation the solution moves on an invariant torus,
% which continuation cannot follow. We can show it instead by the local maxima
% of u from the direct simulations of c5 (Figure 6), which fill an interval.
% Uncomment lines 46-48 to add this information to the plot if you want to.
taurange = [min(A.orbits.tau) - 0.005, max(A.orbits.tau)];
y = [A.hopf_point.u; A.orbits.maxima(:)];
have_sim = isfile(wc_file('data', 'route_to_chaos.mat'));
if have_sim
    C = load(wc_file('data', 'route_to_chaos.mat'), 'diagram');
    %sel = C.diagram.points_tau > A.torus_point.tau & C.diagram.points_tau <= taurange(2);
    %plot(ax, C.diagram.points_tau(sel), C.diagram.points_umax(sel), '.', 'Color', S.torus, 'MarkerSize', 3);
    %y = [y; C.diagram.points_umax(sel)'];
end

mark_point(ax, A.hopf_point.tau, A.hopf_point.u, S.hopf, 'o');
mark_point(ax, A.torus_point.tau, A.torus_point.umax, S.torus, 'd');
xlim(ax, taurange);
% leave a gap below the equilibrium, so it does not sit on the axis
lo = min(y);  hi = max(y);
ylim(ax, [lo - 0.08*(hi - lo), hi + 0.04*(hi - lo)]);
xlabel(ax, '$\tau$');  ylabel(ax, 'maximum of $u$'); % Or local maxima if also plotting data beyond the torus
finish_axes(ax, '(a)');

%% (b) the whole isola of asymmetric orbits and the period-doubled orbits
ax = axes(fig, 'Position', pos{2});  hold(ax, 'on');
draw_isola(ax, B, S);
mark_point(ax, B.isola_events.fold_tau, B.isola_events.fold_umax, S.fold, '^');
mark_point(ax, B.pd_point.tau, B.pd_point.umax, S.pd, 's');
mark_point(ax, B.pd_low.tau, B.pd_low.umax, S.pd, 's');
mark_events(ax, B.ext_events, S);
if ~isempty(B.low_events), mark_events(ax, B.low_events, S); end
all_tau = B.isola_summary.tau;
xlim(ax, [min(all_tau), max(all_tau)] + 0.01*[-1 1]);
xlabel(ax, '$\tau$');  ylabel(ax, 'maximum of $u$');
[h, labels] = legend_items(ax, {'eq', 'po', 'po2', 'hopf', 'torus', 'pd', 'fold'});
if have_sim
    h(end+1) = plot(ax, nan, nan, '.', 'Color', S.torus, 'MarkerSize', 8);
    labels{end+1} = 'torus (simulation)';
end
finish_axes(ax, '(b)');
shade(ax, [0.2655 0.28]);

%% (c) enlargement: the period doubling at tau = 0.276 and the doubled branch
ax = axes(fig, 'Position', pos{3});  hold(ax, 'on');
draw_isola(ax, B, S);
mark_point(ax, B.pd_point.tau, B.pd_point.umax, S.pd, 's');
mark_point(ax, B.folds.tau, B.folds.umax, S.fold, '^');
mark_point(ax, B.pd2.tau, B.pd2.umax, S.pd, 's');
xlim(ax, [0.2655 0.28]);
fit_ylim(ax, B, [0.2655 0.28]);
xlabel(ax, '$\tau$');  ylabel(ax, 'maximum of $u$');
finish_axes(ax, '(c)');

%% (d) enlargement around the folds of the period-doubled orbit
% The window spans the folds and the second period doubling.
events = [B.folds.tau(:); B.pd2.tau(:)];
window = [min(events) - 1.2e-4, max(events) + 1.0e-4];
shade(ax, window);

inside = B.doubled.tau >= window(1) & B.doubled.tau <= window(2);
ax = axes(fig, 'Position', pos{4});  hold(ax, 'on');
plot_branch(ax, B.doubled.tau, B.doubled.umax, B.doubled.stable, S.po2);
mark_point(ax, B.folds.tau, B.folds.umax, S.fold, '^');
mark_point(ax, B.pd2.tau, B.pd2.umax, S.pd, 's');
yz = B.doubled.umax(inside);
xlim(ax, window);
ylim(ax, [min(yz), max(yz)] + 0.1*(max(yz) - min(yz))*[-1 1]);
ax.XAxis.TickValues = linspace(window(1), window(2), 3);
xtickformat(ax, '%.4f');
xlabel(ax, '$\tau$');  ylabel(ax, 'maximum of $u$');
finish_axes(ax, '(d)');

legend_below(fig, h, labels, 5);
save_figure(fig, 'fig2_rich_bifurcations.pdf');

%% ------------------------------------------------------------------------
function draw_isola(ax, B, S)
% The isola (blue, closed up by repeating its first point) and every
% computed branch of period-doubled orbits (orange), each drawn as the
% maximum of u over the orbit, so that every branch is a single curve.
I = B.isola_summary;
close_up = @(x) [x; x(1, :)];
plot_branch(ax, close_up(I.tau), close_up(I.umax), close_up(I.stable), S.po);
branches = {B.doubled, B.doubled_ext, B.doubled_low};
for b = 1:numel(branches)
    D = branches{b};
    if isempty(D), continue, end
    plot_branch(ax, D.tau, D.umax, D.stable, S.po2);
end
end

function mark_events(ax, ev, S)
% Folds and period doublings of a period-doubled branch. Folds closer
% together than 1e-3 in tau would print as one blob, so only the first of
% each cluster is marked; panel (d) shows the cluster near tau = 0.268.
if ~isempty(ev.fold_tau)
    [sorted, order] = sort(ev.fold_tau(:));
    idx = order([true; diff(sorted) > 1e-3]);
    mark_point(ax, ev.fold_tau(idx), ev.fold_umax(idx), S.fold, '^');
end
mark_point(ax, ev.pd_tau, ev.pd_umax, S.pd, 's');
end

function fit_ylim(ax, B, xw)
% y-limits fitted to the curves (maximum of u) inside the window xw.
y = [];
D = {B.isola_summary, B.doubled, B.doubled_ext, B.doubled_low};
for b = 1:numel(D)
    if isempty(D{b}), continue, end
    y = [y; D{b}.umax(D{b}.tau >= xw(1) & D{b}.tau <= xw(2))]; %#ok<AGROW>
end
ylim(ax, [min(y), max(y)] + 0.05*(max(y) - min(y))*[-1 1]);
end

function shade(ax, window)
% Shade the window of the next enlargement, behind everything else.
yl = ylim(ax);
band = patch(ax, window([1 2 2 1]), yl([1 1 2 2]), -ones(1, 4), ...
             [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(band, 'bottom');
ax.SortMethod = 'childorder';
ylim(ax, yl);
end