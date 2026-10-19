%% WC_SETUP  Put this project and DDE-BIFTOOL on the MATLAB path.
%
% Run this once at the start of each MATLAB session, from this folder:
%
%     >> wc_setup
%
% EDIT THE NEXT LINE so that it points to your copy of DDE-BIFTOOL: the
% folder that contains the subfolders 'ddebiftool', 'ddebiftool_extra_psol'
% and 'ddebiftool_utilities'. DDE-BIFTOOL is free to download from
% https://sourceforge.net/projects/ddebiftool/ (or its GitHub mirror).
% The figure scripts only need DDE-BIFTOOL if you rerun the computations.

ddebiftool_folder = fullfile('~', 'ddebiftool');

project = fileparts(mfilename('fullpath'));
addpath(fullfile(project, 'model'), fullfile(project, 'tools'), ...
        fullfile(project, 'compute'), fullfile(project, 'figures'), ...
        fullfile(project, 'listings'));

subfolders = {'ddebiftool', 'ddebiftool_extra_psol', 'ddebiftool_utilities'};
if isfolder(fullfile(ddebiftool_folder, 'ddebiftool'))
    for k = 1:numel(subfolders)
        addpath(fullfile(ddebiftool_folder, subfolders{k}));
    end
    fprintf('Project and DDE-BIFTOOL added to the path.\n');
else
    fprintf(['Project added to the path. DDE-BIFTOOL was not found in\n  %s\n', ...
             'Edit ddebiftool_folder in wc_setup.m to rerun the computations\n', ...
             '(the figure scripts work without it).\n'], ddebiftool_folder);
end

for folder = {'data', 'pdf'}
    if ~isfolder(fullfile(project, folder{1}))
        mkdir(fullfile(project, folder{1}));
    end
end
clear project subfolders k folder ddebiftool_folder
