function plot_branch(ax, x, y, stable, colour)
%PLOT_BRANCH  Draw a branch: solid where stable, dashed where unstable.
%
%   plot_branch(ax, x, y, stable, colour)
%
%   x, y     coordinates of the points along the branch
%   stable   logical array, true where the solution is stable
%   colour   RGB colour (see paper_style)
%
%   The branch is split into stretches of constant stability. Each stretch
%   is extended by one point so that consecutive stretches join up.

S = paper_style();
x = x(:);  y = y(:);  stable = logical(stable(:));
n = numel(x);
i = 1;
while i <= n
    j = i;
    while j < n && stable(j+1) == stable(i)
        j = j + 1;
    end
    last = min(j + 1, n);
    if stable(i), style = '-'; else, style = '--'; end
    plot(ax, x(i:last), y(i:last), style, 'Color', colour, 'LineWidth', S.lw);
    i = j + 1;
end
end
