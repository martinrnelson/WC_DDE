function finish_axes(ax, label)
%FINISH_AXES  Common finishing touches for one panel of a figure.
%
%   finish_axes(ax)          box, outward ticks, thin axis lines
%   finish_axes(ax, label)   also put the panel label, e.g. '(a)', on top

box(ax, 'on');
set(ax, 'TickDir', 'out', 'LineWidth', 0.6, 'Layer', 'top');
if nargin > 1
    title(ax, label, 'FontWeight', 'normal');
end
end
