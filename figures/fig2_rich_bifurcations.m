%% FIGURE 2  Secondary bifurcations of periodic orbits, rich parameter set
%
% (a) Hopf bifurcation of the equilibrium and torus bifurcation of the
%     anti-phase orbits
% (b) period doubling of the asymmetric orbit (blue) and the period-doubled
%     branch (orange) with its folds and a further period doubling
% (c) enlargement of the shaded region in (b)
% Curves show the local maxima of u over one period.
% Data: data/rich_hopf_torus.mat (c2) and data/rich_period_doubling.mat (c3)

clear; close all
A = load(wc_file('data', 'rich_hopf_torus.mat'));
B = load(wc_file('data', 'rich_period_doubling.mat'));
S = paper_style();

fig = new_figure(1.25*S.width, 6);

%% (a) equilibrium, anti-phase orbits, Hopf and torus points
ax = subplot(1, 3, 1);  hold(ax, 'on');
plot_branch(ax, A.equilibrium.tau, A.equilibrium.u, A.equilibrium.stable, S.eq);
for k = 1:size(A.orbits.maxima, 2)        % one maximum per period here
    plot_branch(ax, A.orbits.tau, A.orbits.maxima(:, k), A.orbits.stable, S.po);
end
mark_point(ax, A.hopf_point.tau, A.hopf_point.u, S.hopf, 'o');
mark_point(ax, A.torus_point.tau, A.torus_point.umax, S.torus, 'd');
xlim(ax, [min(A.orbits.tau) - 0.005, max(A.orbits.tau)]);
xlabel(ax, '$\tau$');  ylabel(ax, 'local maxima of $u$');
finish_axes(ax, '(a)');

%% (b) period doubling, period-doubled branch and its folds
ax = subplot(1, 3, 2);  hold(ax, 'on');
for k = 1:size(B.period_T.maxima, 2)      % two maxima per period
    plot_branch(ax, B.period_T.tau, B.period_T.maxima(:, k), B.period_T.stable, S.po);
end
for k = 1:size(B.doubled.maxima, 2)       % four maxima per (doubled) period
    plot_branch(ax, B.doubled.tau, B.doubled.maxima(:, k), B.doubled.stable, S.po2);
end
mark_point(ax, B.pd_point.tau, B.pd_point.maxima, S.pd, 's');
mark_point(ax, B.folds.tau, B.folds.maxima, S.fold, '^');
mark_point(ax, B.pd2.tau, B.pd2.maxima, S.pd, 's');
% the unstable branches stop where the computation stopped (tau = 0.265),
% so start the axis just after that and let them run off the edge
xlim(ax, [0.2655 0.28]);
xlabel(ax, '$\tau$');  ylabel(ax, 'local maxima of $u$');
[h, labels] = legend_items(ax, {'eq', 'po', 'po2', 'hopf', 'torus', 'pd', 'fold'});
finish_axes(ax, '(b)');

%% (c) enlargement around the folds
% The window spans the folds and the second period doubling. Of the four
% maxima, we show the one that changes most inside the window.
events = [B.folds.tau(:); B.pd2.tau(:)];
window = [min(events) - 1.2e-4, max(events) + 1.0e-4];

% shade the window in panel (b), behind everything else
yl = ylim(ax);
band = patch(ax, window([1 2 2 1]), yl([1 1 2 2]), -ones(1, 4), ...
             [0.93 0.93 0.93], 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(band, 'bottom');
ax.SortMethod = 'childorder';
ylim(ax, yl);

inside = B.doubled.tau >= window(1) & B.doubled.tau <= window(2);
change = max(B.doubled.maxima(inside, :)) - min(B.doubled.maxima(inside, :));
[~, kz] = max(change);

ax = subplot(1, 3, 3);  hold(ax, 'on');
plot_branch(ax, B.doubled.tau, B.doubled.maxima(:, kz), B.doubled.stable, S.po2);
mark_point(ax, B.folds.tau, B.folds.maxima(:, kz), S.fold, '^');
mark_point(ax, B.pd2.tau, B.pd2.maxima(:, kz), S.pd, 's');
yz = B.doubled.maxima(inside, kz);
xlim(ax, window);
ylim(ax, [min(yz), max(yz)] + 0.1*(max(yz) - min(yz))*[-1 1]);
ax.XAxis.TickValues = linspace(window(1), window(2), 3);
xtickformat(ax, '%.4f');
xlabel(ax, '$\tau$');  ylabel(ax, 'local maximum of $u$');
finish_axes(ax, '(c)');

legend_below(fig, h, labels, 5);
save_figure(fig, 'fig2_rich_bifurcations.pdf');
