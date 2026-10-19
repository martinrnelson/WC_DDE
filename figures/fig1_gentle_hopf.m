%% FIGURE 1  Delay-induced Hopf bifurcation, gentle parameter set
%
% (a) the equilibrium (black) and the maximum and minimum of u on the
%     periodic orbits (blue), against tau
% (b) the period of the orbits
% Data: data/gentle_hopf.mat, made by compute/c1_gentle_hopf.m

clear; close all
load(wc_file('data', 'gentle_hopf.mat'), 'equilibrium', 'orbits', 'hopf_point');
S = paper_style();

fig = new_figure(S.width, 6);

% (a) bifurcation diagram
ax = subplot(1, 2, 1);  hold(ax, 'on');
plot_branch(ax, equilibrium.tau, equilibrium.u, equilibrium.stable, S.eq);
plot_branch(ax, orbits.tau, orbits.umax, orbits.stable, S.po);
plot_branch(ax, orbits.tau, orbits.umin, orbits.stable, S.po);
mark_point(ax, hopf_point.tau, hopf_point.u, S.hopf, 'o');
xlim(ax, [0 max(orbits.tau)]);
xlabel(ax, '$\tau$');  ylabel(ax, '$u$');
[h, labels] = legend_items(ax, {'eq', 'po', 'hopf'});
finish_axes(ax, '(a)');

% (b) period: the orbit is born with period 2*pi/omega at the Hopf point
ax = subplot(1, 2, 2);  hold(ax, 'on');
plot_branch(ax, orbits.tau, orbits.period, orbits.stable, S.po);
mark_point(ax, hopf_point.tau, 2*pi/hopf_point.omega, S.hopf, 'o');
xlabel(ax, '$\tau$');  ylabel(ax, 'period $T$');
finish_axes(ax, '(b)');

legend_below(fig, h, labels, 5);
save_figure(fig, 'fig1_gentle_hopf.pdf');
