%% FIGURE 3  Floquet multipliers just past the torus and period-doubling points
%
% (a) anti-phase orbit just past the torus bifurcation: a complex pair
%     has left the unit circle
% (b) asymmetric orbit just past the period doubling: a real multiplier
%     has passed through -1
% Filled symbols lie outside the unit circle; the open symbol at 1 is the
% trivial multiplier.
% Data: data/rich_hopf_torus.mat (c2) and data/rich_period_doubling.mat (c3)

clear; close all
A = load(wc_file('data', 'rich_hopf_torus.mat'), 'torus_point');
B = load(wc_file('data', 'rich_period_doubling.mat'), 'pd_point');
S = paper_style();

fig = new_figure(0.75*S.width, 6);

ax = subplot(1, 2, 1);  hold(ax, 'on');
plot_multipliers(ax, A.torus_point.mu_after, S.torus);
finish_axes(ax, sprintf('(a) $\\tau = %.4f$', A.torus_point.tau_after));

ax = subplot(1, 2, 2);  hold(ax, 'on');
plot_multipliers(ax, B.pd_point.mu_after, S.pd);
finish_axes(ax, sprintf('(b) $\\tau = %.4f$', B.pd_point.tau_after));

save_figure(fig, 'fig3_floquet_multipliers.pdf');

%% ------------------------------------------------------------------------
function plot_multipliers(ax, mu, colour)
% Multipliers in the complex plane, with the unit circle.
S = paper_style();
theta = linspace(0, 2*pi, 400);
plot(ax, cos(theta), sin(theta), 'k-', 'LineWidth', 0.6);
outside = abs(mu) > 1 + 1e-6;
plot(ax, real(mu(~outside)), imag(mu(~outside)), 'o', ...
     'Color', [0.4 0.4 0.4], 'MarkerSize', S.ms - 1);
plot(ax, real(mu(outside)), imag(mu(outside)), 'o', ...
     'MarkerFaceColor', colour, 'MarkerEdgeColor', 'k', 'MarkerSize', S.ms);
axis(ax, 'equal');
r = max(1.2, 1.1*max(abs(mu)));
xlim(ax, [-r r]);  ylim(ax, [-r r]);
xlabel(ax, '$\mathrm{Re}\,\mu$');  ylabel(ax, '$\mathrm{Im}\,\mu$');
end
