%% FIGURE 7  Solutions along the route to chaos, rich parameter set
%
% One column per solution: periodic, quasi-periodic, mode-locked, chaotic.
% Top:    u(t) (blue) and v(t) (grey) over a short time window
% Middle: the solution in the (u, v) plane
% Bottom: return map of successive maxima of u. A periodic solution gives
%         a single point, a quasi-periodic one a closed curve, a
%         mode-locked one a finite set of points, a chaotic one a folded set.
% Data: data/route_to_chaos.mat, made by compute/c5_route_to_chaos.m

clear; close all
load(wc_file('data', 'route_to_chaos.mat'), 'examples');
S = paper_style();
n = numel(examples.tau);

fig = new_figure(17, 12);
for k = 1:n
    ts = examples.series(k);
    ax = subplot(3, n, k);  hold(ax, 'on');
    plot(ax, ts.t, ts.u, '-', 'Color', S.po, 'LineWidth', 1);
    plot(ax, ts.t, ts.v, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    xlim(ax, [0 max(ts.t)]);  xlabel(ax, '$t$');
    if k == 1
        ylabel(ax, 'activity');
        h = [plot(ax, nan, nan, '-', 'Color', S.po, 'LineWidth', 1), ...
             plot(ax, nan, nan, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1)];
    end
    finish_axes(ax, {sprintf('(%s) %s', char('a' + k - 1), examples.label{k}), ...
                     sprintf('$\\tau = %.4g$', examples.tau(k))});

    ph = examples.phase(k);
    ax = subplot(3, n, n + k);
    plot(ax, ph.u, ph.v, '-', 'Color', S.po, 'LineWidth', 0.4);
    axis(ax, 'tight');  xlabel(ax, '$u$');
    if k == 1, ylabel(ax, '$v$'); end
    finish_axes(ax);

    m = examples.maxima{k};
    ax = subplot(3, n, 2*n + k);
    plot(ax, m(1:end-1), m(2:end), 'k.', 'MarkerSize', 3);
    lo = min(m);  hi = max(m);  pad = max(0.05*(hi - lo), 1e-4);
    axis(ax, 'equal');  xlim(ax, [lo-pad, hi+pad]);  ylim(ax, [lo-pad, hi+pad]);
    xlabel(ax, '$u_{\max}(n)$');
    if k == 1, ylabel(ax, '$u_{\max}(n+1)$'); end
    finish_axes(ax);
end
legend_below(fig, h, {'$u$', '$v$'}, 2);
save_figure(fig, 'fig7_chaos_examples.pdf');
