$ErrorActionPreference = 'Stop'

$sourceDir = Join-Path $PSScriptRoot '..\src'
$buildDir = Join-Path $PSScriptRoot 'build'
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
Set-Location -LiteralPath $buildDir

$compileFlags = @(
    '-O0', '-g', '-fcheck=all', '-ffree-line-length-none',
    '-fopenmp'
)
$modules = @(
    'config', 'data1d', 'grid', 'gastheory', 'gravity', 'velocity',
    'flux', 'initialization', 'solver', 'output'
)

function Invoke-Native {
    param([string[]]$Arguments)
    & gfortran @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "gfortran failed with exit code $LASTEXITCODE"
    }
}

foreach ($name in $modules) {
    Invoke-Native ($compileFlags + @(
        '-c', (Join-Path $sourceDir "$name.f90"), '-o', "$name.o"
    ))
}

$tests = @(
    @{
        Name = 'test_gauss_hermite'
        Objects = @('config.o', 'data1d.o', 'velocity.o')
    },
    @{
        Name = 'test_force_moments'
        Objects = @('config.o', 'data1d.o', 'gastheory.o', 'gravity.o', 'velocity.o')
    },
    @{
        Name = 'test_temperature_convention'
        Objects = @('config.o', 'data1d.o', 'gastheory.o')
    },
    @{
        Name = 'test_diffuse_wall'
        Objects = @('config.o', 'data1d.o', 'grid.o', 'gastheory.o', 'gravity.o',
                    'velocity.o', 'flux.o')
    },
    @{
        Name = 'test_specular_wall'
        Objects = @('config.o', 'data1d.o', 'grid.o', 'gastheory.o', 'gravity.o',
                    'velocity.o', 'flux.o')
    },
    @{
        Name = 'test_macro_output'
        Objects = @('config.o', 'data1d.o', 'grid.o', 'gastheory.o', 'gravity.o',
                    'velocity.o', 'output.o')
    },
    @{
        Name = 'test_steady_state'
        Objects = @('config.o', 'data1d.o', 'grid.o', 'gastheory.o', 'gravity.o',
                    'velocity.o', 'flux.o', 'solver.o')
    },
    @{
        Name = 'test_primitive_change'
        Objects = @('config.o', 'data1d.o', 'grid.o', 'gastheory.o', 'gravity.o',
                    'velocity.o', 'flux.o', 'solver.o')
    }
)

foreach ($test in $tests) {
    Invoke-Native ($compileFlags + @(
        '-c', (Join-Path $PSScriptRoot "$($test.Name).f90"),
        '-o', "$($test.Name).o"
    ))
    Invoke-Native (@('-fopenmp') + $test.Objects + @(
        "$($test.Name).o", '-o', "$($test.Name).exe"
    ))
    & ".\$($test.Name).exe"
    if ($LASTEXITCODE -ne 0) {
        throw "$($test.Name) failed with exit code $LASTEXITCODE"
    }
}

Invoke-Native ($compileFlags + @(
    '-c', (Join-Path $sourceDir 'main.f90'), '-o', 'main.o'
))
Invoke-Native (@('-fopenmp') + ($modules | ForEach-Object { "$_.o" }) + @(
    'main.o', '-o', 'ugks1d_steady_smoke.exe'
))
& (Join-Path $PSScriptRoot 'test_fourier_steady_smoke.ps1') `
    -Executable (Join-Path $buildDir 'ugks1d_steady_smoke.exe')

Write-Host 'PASS: all UGKS1D regression tests'
