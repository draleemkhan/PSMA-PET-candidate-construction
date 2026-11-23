param(
    [ValidateSet("PET","CT")]
    [string]$Modality = "PET",

    [ValidateSet("internal","external","all")]
    [string]$Scope = "internal",

    [ValidateRange(1,32)]
    [int]$Workers = 4,

    [switch]$RetryFailures
)

$ErrorActionPreference = "Stop"

$ProjectCode = $PSScriptRoot
$PythonScript = Join-Path $ProjectCode "PSMA_objective9D9E_extract_radiomics_DIRECT_RESUMABLE_LONG_PATH_FIXED.py"

if (-not (Test-Path $PythonScript)) {
    throw "Python script not found: $PythonScript"
}

Set-Location $ProjectCode

Write-Host "============================================================"
Write-Host "PSMA Objective 9D/9E direct radiomics extraction"
Write-Host "Modality : $Modality"
Write-Host "Scope    : $Scope"
Write-Host "Workers  : $Workers"
Write-Host "============================================================"

$arguments = @(
    $PythonScript
    "--modality", $Modality
    "--scope", $Scope
    "--workers", $Workers
)

if ($RetryFailures) {
    $arguments += "--retry-failures"
}

python @arguments

if ($LASTEXITCODE -ne 0) {
    throw "Python extraction failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Extraction process finished."
