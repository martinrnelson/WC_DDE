function save_figure(fig, name)
%SAVE_FIGURE  Save a figure as a vector PDF in the pdf/ folder.
%
%   save_figure(fig, 'fig1_gentle_hopf.pdf')
%
%   Uses exportgraphics, available from MATLAB R2020a. On older versions,
%   replace the call with  print(fig, file, '-dpdf', '-painters').

file = wc_file('pdf', name);
exportgraphics(fig, file, 'ContentType', 'vector');
fprintf('Saved %s\n', file);
end
