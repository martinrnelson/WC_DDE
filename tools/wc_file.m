function f = wc_file(folder, name)
%WC_FILE  Full path of a file in one of this project's folders.
%
%   wc_file('data', 'gentle_hopf.mat')   ->  <project>/data/gentle_hopf.mat
%   wc_file('pdf',  'fig1_gentle_hopf.pdf')
%
%   This lets every script read and write the same files, whichever folder
%   MATLAB is currently in.

here = fileparts(mfilename('fullpath'));    % <project>/tools
root = fileparts(here);                     % <project>
f    = fullfile(root, folder, name);
end
