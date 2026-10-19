function fig = new_figure(width, height)
%NEW_FIGURE  A blank figure of the given size (in cm) with LaTeX labels.
%
%   fig = new_figure(width, height)
%
%   Text, tick labels and legends in this figure are typeset with LaTeX,
%   so labels such as '$\tau$' match the paper. The setting applies to
%   this figure only, so other MATLAB figures are unaffected.

S = paper_style();
fig = figure('Units', 'centimeters', 'Position', [2 2 width height], ...
             'Color', 'w');
set(fig, 'defaultTextInterpreter', 'latex', ...
         'defaultAxesTickLabelInterpreter', 'latex', ...
         'defaultLegendInterpreter', 'latex', ...
         'defaultAxesFontSize', S.fs);
end
