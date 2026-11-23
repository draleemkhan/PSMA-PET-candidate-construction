param(
    [Parameter(Mandatory=$true)]
    [string]$PetParamsYaml,

    [Parameter(Mandatory=$true)]
    [string]$CtParamsYaml,

    [string]$PythonExe = "python",

    [ValidateSet("pet","ct","both")]
    [string]$Mode = "both",

    [int]$SaveEvery = 25,

    [double]$SaveMinutes = 5,

    [switch]$RetryFailures
)

$ErrorActionPreference = "Stop"

$ScriptPath = Join-Path $PSScriptRoot "PSMA_objective8D8E_resumable_pyradiomics_SEPARATE_YAML.py"

if (-not (Test-Path $ScriptPath)) {
    throw "Python script not found: $ScriptPath"
}

if (($Mode -eq "pet" -or $Mode -eq "both") -and -not (Test-Path $PetParamsYaml)) {
    throw "PET PyRadiomics YAML not found: $PetParamsYaml"
}

if (($Mode -eq "ct" -or $Mode -eq "both") -and -not (Test-Path $CtParamsYaml)) {
    throw "CT PyRadiomics YAML not found: $CtParamsYaml"
}

$Arguments = @(
    $ScriptPath,
    "--mode", $Mode,
    "--save-every", $SaveEvery,
    "--save-minutes", $SaveMinutes
)

if ($Mode -eq "pet" -or $Mode -eq "both") {
    $Arguments += @("--pet-params-yaml", $PetParamsYaml)
}

if ($Mode -eq "ct" -or $Mode -eq "both") {
    $Arguments += @("--ct-params-yaml", $CtParamsYaml)
}

if ($RetryFailures) {
    $Arguments += "--retry-failures"
}

Write-Host ""
Write-Host "PSMA Objective 8D/8E resumable PyRadiomics extraction"
Write-Host "Mode: $Mode"
Write-Host "PET parameters: $PetParamsYaml"
Write-Host "CT parameters:  $CtParamsYaml"
Write-Host "Progress is checkpointed and resumes automatically."
Write-Host ""

& $PythonExe @Arguments

if ($LASTEXITCODE -ne 0) {
    throw "Python exited with code $LASTEXITCODE"
}
