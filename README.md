# Hermite-UGKS

[![arXiv: 2507.10021](https://img.shields.io/badge/arXiv-2507.10021-b31b1b.svg)](https://arxiv.org/abs/2507.10021)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Fortran reference implementation accompanying the paper **“An efficient solution algorithm for force-driven continuum and rarefied flows”** by Shuangqing Liu, Zuoxu Li, Yonghao Zhang, and Tianbai Xiao ([arXiv:2507.10021](https://arxiv.org/abs/2507.10021)).

The code combines a finite-volume unified gas-kinetic scheme (UGKS) with a Hermite spectral representation of the force term in the Boltzmann–BGK equation. It includes the one- and two-dimensional numerical examples used to study force-driven flows from continuum to rarefied regimes.

## Features

- BGK-based UGKS flux with force-modified particle trajectories
- Gauss–Hermite velocity quadrature
- Hermite spectral evaluation of the velocity-space force term
- Optional finite-difference force operator for comparison
- Well-balanced treatment of hydrostatic equilibrium
- Specular, diffuse-wall, symmetry, and periodic boundary conditions
- Runtime configuration through Fortran namelist files
- Tecplot-compatible output and Python plotting scripts

## Repository layout

```text
Hermite-UGKS/
├── 1d/                     # One-dimensional solver
│   ├── src/                # Fortran source files
│   ├── cases/              # Hydrostatic, Sod, and Fourier inputs
│   ├── tests/              # 1D regression tests
│   ├── input.namelist      # Default 1D configuration
│   ├── Makefile
│   └── README.md           # Detailed 1D documentation
└── 2d/                     # Two-dimensional solver
    ├── src/                # Fortran source files
    ├── cases/              # Poiseuille and Rayleigh–Taylor inputs
    ├── tests/              # 2D regression tests
    ├── input.namelist      # Default 2D configuration
    ├── Makefile
    └── README.md           # Detailed 2D documentation
```

The 1D and 2D programs are independent executables and must be built from their respective directories.

## Requirements

- GNU Fortran (`gfortran`) or Intel Fortran Compiler Classic (`ifort`)
- OpenMP support
- GNU Make for the Unix-like build commands below
- Optional: Python 3 with NumPy and Matplotlib for plotting
- Optional: PowerShell 7 and `gfortran` for the supplied regression-test scripts; these scripts currently invoke `gfortran` directly

## Quick start

Clone or download the repository, then build the solver you want to use.

### 1D solver

```bash
# From the repository root:
cd 1d
make
./ugks1d cases/sod_kn001_gh34.namelist
```

### 2D solver

```bash
# From the repository root:
cd 2d
make
./ugks2d cases/poiseuille_fast.namelist
```

Run each executable from its solver directory. Without an argument, it reads `input.namelist` from the current working directory:

```bash
./ugks1d
# or
./ugks2d
```

The Makefiles use `gfortran` by default. To build with Intel Fortran instead, clean any files produced by the other compiler and set `FC`:

```bash
make clean
make FC=ifort
```

The Makefiles automatically select the appropriate module-output and OpenMP options for `gfortran` or `ifort`. Always run `make clean` when switching compilers because their object and module files are not interchangeable.

Compiler-specific optimization flags can be supplied explicitly:

```bash
# GNU Fortran
make FFLAGS="-O3 -ffree-line-length-none -fno-range-check"

# Intel Fortran
make FC=ifort FFLAGS="-O3"
```

On Windows with PowerShell, the 2D solver also provides a native build script:

```powershell
cd 2d
./build.ps1
./ugks2d.exe ./cases/poiseuille_fast.namelist
```

## Included benchmarks

| Dimension | Case | `case_id` | Ready-to-run input |
|---|---|---:|---|
| 1D | Hydrostatic equilibrium | 2 | `1d/cases/hydrostatic_kn1_gh36.namelist` |
| 1D | Sod shock tube under an external force | 3 | `1d/cases/sod_kn001_gh34.namelist` |
| 1D | Force-driven Fourier flow | 4 | `1d/cases/fourier_force2_kn01_gh22.namelist` |
| 2D | Force-driven Poiseuille flow | 5 | `2d/cases/poiseuille_fast.namelist` |
| 2D | Rayleigh–Taylor instability | 6 | `2d/cases/rt_60_fast.namelist` or `rt_120.namelist` |

The “fast” 2D inputs are reduced validation runs. The larger or steady configurations can require substantially more time and memory because the distribution function is resolved in both physical and velocity space.

Convenience Make targets are also provided:

```bash
make -C 1d run2          # hydrostatic equilibrium
make -C 1d run3          # Sod shock tube
make -C 1d run4          # Fourier flow
make -C 2d run5          # fast Poiseuille case
make -C 2d run6          # 60 x 60 Rayleigh–Taylor case
```

## Configuration

Each run is controlled by a Fortran `&params` namelist. Important fields include:

| Parameter | Description |
|---|---|
| `nx`, `ny` | Number of cells in physical space (`ny` is 2D only) |
| `xlength`, `ylength` | Physical-domain dimensions |
| `nv` | Number of velocity nodes per velocity dimension |
| `vmax`, `Tref` | Velocity-space extent and reference-temperature scaling |
| `quad_type` | `gauss_hermite` or `newton_cotes` |
| `force_type` | Force-term discretization; see below |
| `Nm` | Hermite expansion order |
| `cfl`, `maxtime` | CFL number and final/safety time |
| `case_id` | Benchmark selector |
| `knudsen` | Reference Knudsen number |
| `phi` | External-force magnitude |
| `TwL`, `TwR` | Wall temperatures for diffuse boundaries |

The 1D solver accepts `finite_difference` and `hermite` force operators. The 2D solver additionally contains `hermite_symmetric` and `hermite_symmetric_b`, two equivalent symmetric-tensor formulations; the default 2D input uses `hermite_symmetric`.

Start from a file in `cases/` when reproducing a benchmark. For example:

```bash
cp cases/hydrostatic_kn1_gh36.namelist my_case.namelist
# edit my_case.namelist, then run from 1d/
./ugks1d my_case.namelist
```

## Output and plotting

The solvers write Tecplot-compatible files named `macro<iteration>.dat` in the current working directory.

- The 1D output contains position, density, velocity, temperature, pressure, heat flux, pressure perturbation, and Hermite diagnostics.
- The 2D output contains the mesh and cell-centered macroscopic fields, including density, velocity, temperature, pressure, and heat flux.

Run simulations from a dedicated directory, or move existing `macro*.dat` files first, to avoid mixing results from different configurations.

Plotting scripts are supplied in both solver directories. They require NumPy and Matplotlib:

```bash
python3 -m pip install numpy matplotlib
python3 plot_hermite_validation.py
```

The scripts expect result files in the directory paths defined near the top of each script. Update those paths to match the output of your run. Figures are written below `figures/`.

## Tests

Regression tests are provided as PowerShell scripts:

```powershell
cd 1d
./tests/run_tests.ps1

cd ../2d
./tests/run_tests.ps1
```

The scripts compile debug builds with runtime checks and exercise velocity quadrature, force moments, boundary conditions, conservation-related behavior, and smoke cases. They create test build artifacts below each `tests/build/` directory.

## Method overview

The governing kinetic equation is the Boltzmann equation with an external acceleration, closed with the BGK relaxation model. The implementation evolves both the particle distribution and its conservative moments. UGKS constructs the interface flux from the integral solution of the BGK equation, coupling free transport and collision over a time step. The acceleration term is evaluated using a Hermite expansion in velocity space, while Gauss–Hermite nodes provide an efficient collocation and quadrature rule for near-Maxwellian distributions.

See Sections 2–4 of the paper for the derivation, numerical analysis, and benchmark definitions. Implementation details specific to each executable are documented in [`1d/README.md`](1d/README.md) and [`2d/README.md`](2d/README.md).

## Citation

If this code is useful in your work, please cite the accompanying paper:

```bibtex
@article{liu2025efficient,
  title   = {An efficient solution algorithm for force-driven continuum and rarefied flows},
  author  = {Liu, Shuangqing and Li, Zuoxu and Zhang, Yonghao and Xiao, Tianbai},
  year    = {2025},
  journal = {arXiv preprint arXiv:2507.10021},
  doi     = {10.48550/arXiv.2507.10021}
}
```

## License

This project is licensed under the [MIT License](LICENSE). You may use, copy, modify, merge, publish, distribute, sublicense, and sell copies of the software, subject to the conditions in the license.
