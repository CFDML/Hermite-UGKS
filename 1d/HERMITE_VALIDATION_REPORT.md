# `ugks1d` Hermite validation against `ugks_hermite.pdf`

Date: 2026-07-28

## Scope

Only the Hermite force discretization with Gauss-Hermite (GH) quadrature was
tested. Newton-Cotes and two-dimensional cases were not tested. In accordance
with the requested scope, one plotted observable was checked for each of the
three one-dimensional cases:

| Case | Paper setting used | Checked observable |
| --- | --- | --- |
| Hydrostatic | `Kn=1`, `Nx=600`, GH36, `Nm=6`, `a=-1`, `t=0.08` | pressure perturbation, Figure 2(d) |
| Forced Sod | `Kn=0.01`, `Nx=600`, GH34, `Nm=6`, `a=-1`, `t=0.2` | density, Figure 6(a) |
| Fourier | `Kn=0.1`, `Nx=800`, GH22, `Nm=10`, `a=2`, steady | density, Figure 7(a) |

The PDF contains plotted curves rather than numerical reference tables.
Consequently, "paper recovery" below means agreement with the values and
shapes readable from those figures, together with independent discrete
moment, boundary-flux, mass, and steady-state tests. It is not claimed as an
independent experimental validation of the physical model.

## Recovery results

| Case and checked value | Corrected GH/Hermite result | Paper figure | Assessment |
| --- | ---: | ---: | --- |
| Hydrostatic, `p-p_ref` near `x=0.5` | `3.31210e-3` | about `3.3e-3` | pass |
| Sod, `rho` near `x=0.4` | `0.638376` | about `0.64` | pass |
| Fourier, left/right `rho` | `0.09797 / 3.72486` | about `0.10 / 3.7` | pass |

The recovered hydrostatic perturbation is
`3.65338e-3, 3.31210e-3, 3.52749e-3` near
`x=0.4, 0.5, 0.6`, reproducing the central dip in Figure 2(d). The Sod density
profile reproduces the Figure 6(a) shock/contact structure. The Fourier
density and temperature profiles reproduce Figure 7(a,b); the final mass is
`1.0000000000000007`.

The Fourier run stopped only after 1000 consecutive accepted steps:

```text
iteration = 213439
time = 18.128677
max fixed-scale primitive change = 7.6870665e-10
relative mass drift = 1.6875390e-14
```

## Errors found and corrected

### Mathematics and numerical discretization

1. **Invalid GH generator.** Independent Newton iterations converged to
   duplicate outer roots and produced an unordered grid. The generator now
   uses a normalized Hermite recurrence, symmetry-paired ordered roots, and
   the standard GH weights. Tests cover ordering, uniqueness, symmetry,
   Maxwellian moments, and force moments.
2. **Temperature convention.** The paper uses `R=1/2`; therefore
   `T=1/lambda` and `p=rho*T/2`. Fourier initialization, wall states, and
   output previously mixed this with `lambda=1/(2T)`. They now use one
   consistent convention.
3. **Steady residual.** The old momentum normalization could divide by a
   vanishing velocity. The new residual uses fixed physical scales for
   density, velocity, and temperature and requires both residual and global
   mass criteria for consecutive steps.

### Boundary physics

1. **Diffuse wall.** The old outgoing denominator did not use the same
   reduced Maxwellian as the boundary distribution and therefore did not
   enforce impermeability. The code now constructs a unit-density wall
   Maxwellian and solves its discrete outgoing coefficient directly:
   `rho_w=-R_in/R_out`. Equilibrium left and right wall mass fluxes are zero
   to roundoff in the regression test.
2. **Specular wall.** The old implementation reflected the entire
   distribution. It now replaces only the outgoing half range and preserves
   the incoming half. An asymmetric-distribution regression test verifies
   zero normal mass flux at both walls.

### Code and execution

1. The Windows Makefile used Unix-only commands; it now builds and cleans
   serially on Windows.
2. The configuration path buffer truncated long case paths; it is now length
   260.
3. Hermite diagnostic variables are initialized, macro output includes
   physical heat flux, and Fourier reports a real steady-state decision.
4. Seven unit/regression programs plus a consecutive-steady-state smoke test
   are provided in `tests/run_tests.ps1`.

## Heat-flux convention found during the audit

Paper Eq.(7) defines

```text
q = integral 0.5*(v-V)*(v-V)^2*f dv
```

and the corrected unified code follows that definition. The supplied legacy
file `算例4.Fourier/fourier_flow_f_1.f90` computes the same reduced heat flux
without the factor `0.5`. Its output is therefore twice the physical heat
flux defined by the paper. This legacy inconsistency was not copied into the
corrected solver and no empirical scale factor is used in the validation
plot.

## Reproduction

Build and run all regression tests:

```powershell
cd ugks1d
make clean
make
.\tests\run_tests.ps1
```

The three run inputs are in `cases/`. After running them, generate the simple
three-panel plot with:

```powershell
C:\Users\ps\anaconda3\python.exe .\plot_hermite_validation.py
```

This writes `figures/hermite_1d_validation.png` and
`figures/hermite_1d_validation.pdf`.
