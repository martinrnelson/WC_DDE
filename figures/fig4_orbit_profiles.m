%% FIGURE 4  Three periodic orbits of the rich parameter set
%
% Top row:    u(t) (coloured) and v(t) (grey) over two periods of the
%             period-T orbit, with the equilibrium value dotted
% Bottom row: the same orbits in the (u, v) plane, with the symmetry line
%             u = v (dotted) and the equilibrium (dot)
% (a) anti-phase orbit:  v(t) = u(t + T/2), so the picture is symmetric
% (b) asymmetric orbit:  two different maxima, no symmetry
% (c) period-doubled orbit: maxima alternate, period 2T
% Data: data/rich_hopf_torus.mat (c2) and data/rich_period_doubling.mat (c3)

clear; close all
A = load(wc_file('data', 'rich_hopf_torus.mat'), 'profile_antiphase', 'equilibrium');
B = load(wc_file('data', 'rich_period_doubling.mat'), 'profile_asym', 'profile_doubled');
S = paper_style();
ustar = A.equilibrium.u(1);           % u* = v* for this parameter set

orbits  = {A.profile_antiphase, B.profile_asym, B.profile_doubled};
titles  = {'anti-phase', 'asymmetric, period $T$', 'period-doubled, period $2T$'};
colours = {S.po, S.po, S.po2};
repeats = [2, 2, 1];   % copies of each orbit needed to fill two periods T

fig = new_figure(S.width, 9);
for k = 1:3
    pr = orbits{k};
    % time in units of T, the period of the period-T orbit: 0 to 2
    s = [];  U = [];  V = [];
    for r = 0:repeats(k)-1
        s = [s, pr.mesh + r];  U = [U, pr.u];  V = [V, pr.v];   %#ok<AGROW>
    end
    s = s * 2/repeats(k);
    T = pr.period * repeats(k)/2;

    % top row: time series
    ax = subplot(2, 3, k);  hold(ax, 'on');
    plot(ax, s, U, '-', 'Color', colours{k}, 'LineWidth', S.lw);
    plot(ax, s, V, '-', 'Color', S.grey, 'LineWidth', S.lw);
    yline(ax, ustar, ':', 'Color', 'k', 'LineWidth', 0.6);
    xlim(ax, [0 2]);  xticks(ax, 0:0.5:2);  xlabel(ax, '$t/T$');
    if k == 1
        ylabel(ax, 'activity');
        h = [plot(ax, nan, nan, '-', 'Color', S.po,   'LineWidth', S.lw), ...
             plot(ax, nan, nan, '-', 'Color', S.grey, 'LineWidth', S.lw)];
    end
    finish_axes(ax, {sprintf('(%s) %s', char('a' + k - 1), titles{k}), ...
                     sprintf('$\\tau = %.3f$, $T = %.3f$', pr.tau, T)});

    % bottom row: (u, v) plane
    ax = subplot(2, 3, 3 + k);  hold(ax, 'on');
    lo = min([U V]);  hi = max([U V]);  pad = 0.05*(hi - lo);
    plot(ax, [lo-pad, hi+pad], [lo-pad, hi+pad], ':', 'Color', S.grey);
    plot(ax, pr.u, pr.v, '-', 'Color', colours{k}, 'LineWidth', S.lw);
    plot(ax, ustar, ustar, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', S.ms - 2);
    axis(ax, 'equal');  xlim(ax, [lo-pad, hi+pad]);  ylim(ax, [lo-pad, hi+pad]);
    % three rounded ticks on each axis, so the labels do not overlap
    decimals = max(0, 1 - floor(log10((hi - lo)/2)));
    ticks = unique(round(linspace(lo, hi, 3), decimals));
    ax.XAxis.TickValues = ticks;  ax.YAxis.TickValues = ticks;
    xtickformat(ax, sprintf('%%.%df', decimals));
    ytickformat(ax, sprintf('%%.%df', decimals));
    xlabel(ax, '$u$');
    if k == 1, ylabel(ax, '$v$'); end
    finish_axes(ax);
end
legend_below(fig, h, {'$u$', '$v$'}, 2);
save_figure(fig, 'fig4_orbit_profiles.pdf');
