function mark_point(ax, x, y, colour, marker)
%MARK_POINT  Mark bifurcation points with a filled, black-edged marker.
%
%   mark_point(ax, x, y, colour, marker)
%
%   Markers used in the paper:  'o' Hopf,  'd' torus,  's' period doubling,
%   '^' fold of periodic orbits,  'p' torus resonance.

S = paper_style();
plot(ax, x, y, marker, 'MarkerFaceColor', colour, 'MarkerEdgeColor', 'k', ...
     'MarkerSize', S.ms, 'LineWidth', 0.5, 'LineStyle', 'none');
end
