$ErrorActionPreference = "Stop"

$CodeDir = $PSScriptRoot
$Runner  = Join-Path $CodeDir "Run-PSMA-Objective9D9E-Direct-Radiomics.ps1"

if (-not (Test-Path $Runner)) {
    throw "Runner not found: $Runner"
}

Set-Location $CodeDir

Write-Host "============================================================"
Write-Host "DEEP-PSMA EXTERNAL CT RADIOMICS"
Write-Host "Scope: external"
Write-Host "Workers: 4"
Write-Host "Resumable: YES"
Write-Host "Periodic saving: YES"
Write-Host "Candidate source: corrected PSMA-only Objective 9A maps"
Write-Host "============================================================"

powershell -ExecutionPolicy Bypass -File $Runner `
    -Modality CT `
    -Scope external `
    -Workers 4

if ($LASTEXITCODE -ne 0) {
    throw "External CT radiomics failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "External CT radiomics process finished."
