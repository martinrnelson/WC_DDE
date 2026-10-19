function [h, labels] = legend_items(ax, items)
%LEGEND_ITEMS  Legend entries for the objects shown in a bifurcation diagram.
%
%   [h, labels] = legend_items(ax, items)
%
%   items is a cell array chosen from
%     'eq' 'po' 'po2'                       (curves)
%     'hopf' 'torus' 'pd' 'fold'            (bifurcation points)
%   The entries 'stable' (solid) and 'unstable' (dashed) are always added
%   first, drawn in grey so they are not mistaken for a particular branch.
%
%   The entries are invisible dummy objects in ax (they have NaN
%   coordinates), so the legend can describe styles independently of the
%   data. Pass h and labels to legend_below.

S = paper_style();
h = [plot(ax, nan, nan, '-',  'Color', S.grey, 'LineWidth', S.lw), ...
     plot(ax, nan, nan, '--', 'Color', S.grey, 'LineWidth', S.lw)];
labels = {'stable', 'unstable'};

curve  = @(c) plot(ax, nan, nan, '-', 'Color', c, 'LineWidth', S.lw);
marker = @(c, m) plot(ax, nan, nan, m, 'MarkerFaceColor', c, ...
                      'MarkerEdgeColor', 'k', 'MarkerSize', S.ms);
for k = 1:numel(items)
    switch items{k}
        case 'eq',    h(end+1) = curve(S.eq);           labels{end+1} = 'equilibrium';            %#ok<AGROW>
        case 'po',    h(end+1) = curve(S.po);           labels{end+1} = 'periodic orbit';         %#ok<AGROW>
        case 'po2',   h(end+1) = curve(S.po2);          labels{end+1} = 'period-doubled orbit';   %#ok<AGROW>
        case 'hopf',  h(end+1) = marker(S.hopf, 'o');   labels{end+1} = 'Hopf';                   %#ok<AGROW>
        case 'torus', h(end+1) = marker(S.torus, 'd');  labels{end+1} = 'torus';                  %#ok<AGROW>
        case 'pd',    h(end+1) = marker(S.pd, 's');     labels{end+1} = 'period doubling';        %#ok<AGROW>
        case 'fold',  h(end+1) = marker(S.fold, '^');   labels{end+1} = 'fold of periodic orbits';%#ok<AGROW>
        otherwise,    error('legend_items: unknown item ''%s''', items{k})
    end
end
end
