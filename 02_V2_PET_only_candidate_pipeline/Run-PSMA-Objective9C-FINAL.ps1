$ErrorActionPreference = "Stop"

$CodeDir = $PSScriptRoot
$Script = Join-Path $CodeDir "PSMA_objective9C_extract_routine_features_FINAL.py"

if (-not (Test-Path $Script)) {
    throw "Script not found: $Script"
}

Set-Location $CodeDir

Write-Host "============================================================"
Write-Host "PSMA OBJECTIVE 9C - FINAL ROUTINE FEATURE EXTRACTION"
Write-Host "Automatically uses/creates O: short-path alias"
Write-Host "Processes internal train/validation/test and external"
Write-Host "Resumable checkpointing after every completed case"
Write-Host "============================================================"

python $Script --cohort all

if ($LASTEXITCODE -ne 0) {
    throw "Objective 9C failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Objective 9C completed successfully."
