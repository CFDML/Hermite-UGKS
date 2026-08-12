$ErrorActionPreference = 'Stop'

$sourceDir = Join-Path $PSScriptRoot '..\src'
$buildDir = Join-Path $PSScriptRoot 'build'
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
Push-Location -LiteralPath $buildDir

$flags = @(
    '-O0', '-g', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
    '-ffree-line-length-none', '-fopenmp'
)
$modules = @(
    'config', 'data2d', 'grid', 'gastheory2d', 'gravity',
    'hermite_symmetric2d', 'velocity',
    'boundary', 'flux', 'initialization', 'solver', 'output2d'
)

function Invoke-Fortran {
    param([string[]]$Arguments)
    & gfortran @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "gfortran failed with exit code $LASTEXITCODE"
    }
}

foreach ($name in $modules) {
    Invoke-Fortran ($flags + @(
        '-J.', '-I.', '-c', (Join-Path $sourceDir "$name.f90"), '-o', "$name.o"
    ))
}

$tests = @(
    @{
        Name = 'test_velocity_force'
        Objects = @('config.o','data2d.o','gastheory2d.o','gravity.o','velocity.o')
    },
    @{
        Name = 'test_symmetric_hermite'
        Objects = @('config.o','data2d.o','gastheory2d.o','gravity.o',
                    'hermite_symmetric2d.o','velocity.o')
    },
    @{
        Name = 'test_boundaries'
        Objects = @('config.o','data2d.o','grid.o','gastheory2d.o','gravity.o',
                    'velocity.o','boundary.o','flux.o')
    },
    @{
        Name = 'test_force_and_ghosts'
        Objects = @('config.o','data2d.o','grid.o','gastheory2d.o','gravity.o',
                    'velocity.o','boundary.o','flux.o')
    }
)

foreach ($test in $tests) {
    Invoke-Fortran ($flags + @(
        '-J.', '-I.', '-c', (Join-Path $PSScriptRoot "$($test.Name).f90"),
        '-o', "$($test.Name).o"
    ))
    Invoke-Fortran (@('-fopenmp') + $test.Objects + @(
        "$($test.Name).o", '-o', "$($test.Name).exe"
    ))
    & ".\$($test.Name).exe"
    if ($LASTEXITCODE -ne 0) {
        throw "$($test.Name) failed with exit code $LASTEXITCODE"
    }
}

Write-Host 'PASS: all UGKS2D regression tests'
Pop-Location
