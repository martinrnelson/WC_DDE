function S = paper_style()
%PAPER_STYLE  Colours, line widths and sizes shared by all figures.
%
%   S = paper_style()
%
%   Conventions used in every bifurcation diagram:
%     solid line = stable,  dashed line = unstable
%     colour     = type of object (Okabe-Ito palette, colour-blind safe)

S.eq    = [0 0 0];              % equilibria                     black
S.po    = [0 114 178]/255;      % periodic orbits                blue
S.po2   = [230 159 0]/255;      % period-doubled orbits          orange
S.hopf  = [213 94 0]/255;       % Hopf points and curves         vermillion
S.torus = [204 121 167]/255;    % torus points and curves        reddish purple
S.pd    = [0 158 115]/255;      % period-doubling points/curves  bluish green
S.fold  = [86 180 233]/255;     % folds of periodic orbits       sky blue
S.grey  = [0.55 0.55 0.55];

S.lw = 1.2;                     % line width
S.ms = 5;                       % marker size
S.fs = 9;                       % font size (points)
S.width = 16;                   % width of a full-width figure (cm)
end
