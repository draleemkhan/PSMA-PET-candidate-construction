% set_psma_paths_TEMPLATE.m
% Copy this file to set_psma_paths.m, edit the three folders and run it once
% per MATLAB session before running any MATLAB script of this repository.
% set_psma_paths.m is listed in .gitignore and is never committed.

setenv("PSMA_PROJECT_ROOT",  "C:\path\to\output_root");          % outputs are written here
setenv("PSMA_INTERNAL_DATA", "C:\path\to\PSMA-PET-CT-Lesions_v3"); % imagesTr, labelsTr, psma_metadata.csv
setenv("PSMA_EXTERNAL_DATA", "C:\path\to\DEEP-PSMA_extracted");    % extracted DEEP-PSMA cases

% Make psma_root.m available to all scripts
addpath(fileparts(mfilename("fullpath")));
