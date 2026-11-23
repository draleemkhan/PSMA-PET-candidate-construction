# set_psma_paths_TEMPLATE.ps1
# Copy to set_psma_paths.ps1, edit the three folders and dot-source it in the
# PowerShell session that runs the Python steps:   . .\config\set_psma_paths.ps1
# The values must match those used in MATLAB (config/set_psma_paths.m).

$env:PSMA_PROJECT_ROOT  = "C:\path\to\output_root"
$env:PSMA_INTERNAL_DATA = "C:\path\to\PSMA-PET-CT-Lesions_v3"
$env:PSMA_EXTERNAL_DATA = "C:\path\to\DEEP-PSMA_extracted"
