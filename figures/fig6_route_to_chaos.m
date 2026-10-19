%% FIGURE 6  Route to chaos through torus breakdown, rich parameter set
%
% (a) orbit diagram: local maxima of u after transients, by simulation
% (b) largest Lyapunov exponent
% Dashed lines mark the Hopf point, the torus point and the onset of chaos.
% Data: data/route_to_chaos.mat (c5) and data/rich_hopf_torus.mat (c2)

clear; close all
load(wc_file('data', 'route_to_chaos.mat'), 'diagram', 'lyap');
A = load(wc_file('data', 'rich_hopf_torus.mat'), 'hopf_point', 'torus_point');
S = paper_style();
marks   = [A.hopf_point.tau, A.torus_point.tau, lyap.tau_chaos];
mcolour = {S.hopf, S.torus, [0.35 0.35 0.35]};
mlabel  = {'$\tau_H$', '$\tau_T$', '$\tau_C$'};

fig = new_figure(14, 10);

ax1 = subplot(2, 1, 1);  hold(ax1, 'on');
plot(ax1, diagram.points_tau, diagram.points_umax, 'k.', 'MarkerSize', 1.5);
ylabel(ax1, 'local maxima of $u$');
finish_axes(ax1, '(a)');

ax2 = subplot(2, 1, 2);  hold(ax2, 'on');
yline(ax2, 0, ':', 'Color', [0.5 0.5 0.5]);
plot(ax2, lyap.tau, lyap.lambda, '.-', 'Color', S.po, 'LineWidth', 1, 'MarkerSize', 6);
xlabel(ax2, '$\tau$');  ylabel(ax2, '$\lambda_{\max}$');
finish_axes(ax2, '(b)');

for ax = [ax1 ax2]
    for k = 1:3
        xline(ax, marks(k), '--', 'Color', mcolour{k}, 'LineWidth', 0.8);
    end
    xlim(ax, [diagram.tau(1), diagram.tau(end)]);
end
yl = ylim(ax1);
for k = 1:3
    text(ax1, marks(k), yl(2), ['\,' mlabel{k}], 'Color', mcolour{k}, ...
         'VerticalAlignment', 'top');
end
linkaxes([ax1 ax2], 'x');
save_figure(fig, 'fig6_route_to_chaos.pdf');
