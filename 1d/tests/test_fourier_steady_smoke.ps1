param(
    [Parameter(Mandatory = $true)]
    [string]$Executable
)

$config = Join-Path $PSScriptRoot 'fourier_steady_smoke.namelist'
$output = & $Executable $config 2>&1
if ($LASTEXITCODE -ne 0) {
    $output
    throw "Fourier steady smoke executable failed with exit code $LASTEXITCODE"
}
if (($output -join "`n") -notmatch 'Fourier steady state reached:\s+3') {
    $output
    throw 'Fourier steady state was not held for the requested consecutive steps'
}
$output
Write-Host 'PASS: Fourier main loop uses the steady-state predicate'
