%% FIGURE 5  Two-parameter bifurcation diagram in the (b, tau) plane
%
% Hopf curves of the equilibrium (exact, Section 3.1): anti-phase solid,
% in-phase dash-dot; the equilibrium is stable in the shaded region; dots
% mark double Hopf points. The torus curve (purple) and the period-doubling
% curve (green, a closed loop) come from DDE-BIFTOOL. Stars mark rotation
% numbers 1/2, 1/3 and 1/4 on the torus curve. The dotted line is the
% slice b = 2.5 of Figures 2-7.
% Data: data/two_parameter.mat, made by compute/c4_two_parameter.m

clear; close all
load(wc_file('data', 'two_parameter.mat'), 'p', 'hopf', 'double_hopf', ...
     'torus', 'pd_curve', 'resonance');
S = paper_style();

fig = new_figure(10, 8);
ax = axes(fig);  hold(ax, 'on');

% the equilibrium is stable below the first anti-phase Hopf curve
patch(ax, [hopf.b, fliplr(hopf.b)], [hopf.tau_antiphase(1,:), zeros(size(hopf.b))], ...
      [0.92 0.92 0.92], 'EdgeColor', 'none');
plot(ax, hopf.b, hopf.tau_antiphase, '-',  'Color', S.hopf, 'LineWidth', S.lw);
plot(ax, hopf.b, hopf.tau_inphase,   '-.', 'Color', S.hopf, 'LineWidth', S.lw);
plot(ax, double_hopf.b, double_hopf.tau, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', S.ms - 1);

plot(ax, torus.b,    torus.tau,    '-', 'Color', S.torus, 'LineWidth', S.lw);
plot(ax, pd_curve.b, pd_curve.tau, '-', 'Color', S.pd,    'LineWidth', S.lw);

mark_point(ax, resonance.b, resonance.tau, S.torus, 'p');
for k = 1:numel(resonance.b)
    text(ax, resonance.b(k) + 0.08, resonance.tau(k) - 0.004, ...
         sprintf('$1/%d$', resonance.n(k)), 'VerticalAlignment', 'top', ...
         'FontSize', S.fs - 1, 'BackgroundColor', 'w', 'Margin', 0.5);
end

xline(ax, p.b, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
xlim(ax, [min(hopf.b), max(hopf.b)]);  ylim(ax, [0 0.35]);
xlabel(ax, '$b = c$');  ylabel(ax, '$\tau$');
finish_axes(ax);

% legend entries (invisible dummy objects with the right styles)
h = [plot(ax, nan, nan, '-',  'Color', S.hopf,  'LineWidth', S.lw), ...
     plot(ax, nan, nan, '-.', 'Color', S.hopf,  'LineWidth', S.lw), ...
     plot(ax, nan, nan, '-',  'Color', S.torus, 'LineWidth', S.lw), ...
     plot(ax, nan, nan, '-',  'Color', S.pd,    'LineWidth', S.lw), ...
     plot(ax, nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', S.ms - 1), ...
     plot(ax, nan, nan, 'p', 'MarkerFaceColor', S.torus, 'MarkerEdgeColor', 'k', 'MarkerSize', S.ms), ...
     patch(ax, nan, nan, [0.92 0.92 0.92], 'EdgeColor', 'none')];
labels = {'Hopf (anti-phase)', 'Hopf (in-phase)', 'torus', 'period doubling', ...
          'double Hopf', 'torus resonance', 'stable equilibrium'};
legend_below(fig, h, labels, 3);
save_figure(fig, 'fig5_two_parameter.pdf');
