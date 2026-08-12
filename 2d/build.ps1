$ErrorActionPreference = 'Stop'

$sourceDir = Join-Path $PSScriptRoot 'src'
$buildDir = Join-Path $PSScriptRoot 'build'
$executable = Join-Path $PSScriptRoot 'ugks2d.exe'
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
Push-Location -LiteralPath $buildDir

try {
    $flags = @('-O2','-ffree-line-length-none','-fopenmp')
    $modules = @(
        'config','data2d','grid','gastheory2d','gravity','hermite_symmetric2d','velocity',
        'boundary','flux','initialization','solver','output2d','main'
    )

    foreach ($name in $modules) {
        & gfortran @flags '-J.' '-I.' '-c' (Join-Path $sourceDir "$name.f90") '-o' "$name.o"
        if ($LASTEXITCODE -ne 0) {
            throw "Compilation failed for $name.f90"
        }
    }

    $objects = $modules | ForEach-Object { "$_.o" }
    & gfortran '-fopenmp' @objects '-o' $executable
    if ($LASTEXITCODE -ne 0) {
        throw 'Link failed'
    }
    Write-Host "Built $executable"
}
finally {
    Pop-Location
}
