function legend_below(fig, h, labels, ncol)
%LEGEND_BELOW  One legend for the whole figure, in a strip below the panels.
%
%   legend_below(fig, h, labels, ncol)
%
%   The figure is made taller by the height of the strip and every panel
%   is moved up by the same amount, so the legend never covers any data.

S = paper_style();
panels = findobj(fig, 'Type', 'axes');
set(panels, 'Units', 'centimeters');
nrows  = ceil(numel(h) / ncol);
strip  = 0.45*nrows + 0.35;                        % strip height (cm)

fig.Units = 'centimeters';
pos = fig.Position;  pos(4) = pos(4) + strip;  fig.Position = pos;
for ax = panels(:)'
    p = ax.Position;  p(2) = p(2) + strip;  ax.Position = p;
end

lgd = legend(ancestor(h(1), 'axes'), h, labels, 'Orientation', 'horizontal', ...
             'NumColumns', ncol, 'Box', 'off', 'FontSize', S.fs, ...
             'Units', 'centimeters', 'AutoUpdate', 'off');
lp = lgd.Position;
lgd.Position = [(pos(3) - lp(3))/2, 0.1, lp(3), lp(4)];
end
