# UGKS2D Hermite validation

This note records what was checked for the two requested two-dimensional
examples.  It does not validate the Newton--Cotes branch.

## Reference settings

| case | physical grid | velocity rule | force expansion | other parameters |
|---|---:|---:|---:|---|
| Poiseuille | 3 x 100 | GH16 | Hermite, Nm=7 | Kn=0.02, ax=0.1, ay=0, Pr=1 |
| Rayleigh--Taylor fast | 60 x 60 | GH20 | Hermite, Nm=1 | Kn=0.01, radial phi=-1.5, Pr=2/3 |
| Rayleigh--Taylor check | 120 x 120 | GH20 | Hermite, Nm=1 | Kn=0.01, radial phi=-1.5, Pr=2/3 |

The RT domain is the square quarter-domain `[0,1] x [0,1]`.  The left and
bottom boundaries are symmetry axes.  The right and top boundaries are closed
specular faces, matching the mass-conserving effect of the legacy code (which
does not compute an outer-face flux and instead replaces its last physical
row/column).  A literal zero-gradient outer-face implementation was rejected:
it changed the total mass by `3.46e-5` in only `t=0.001`, whereas the specular
implementation changed it by `4.60e-16` in the same debug smoke test.

## Source audit and corrections

The unified solver was compared line by line with `算例5poiseuille` and
`算例6.Rayleigh-Taylor`.  The main corrections are:

1. Generate distinct, ordered Gauss--Hermite abscissae using the normalized
   Hermite-root recurrence and convert the weighted GH weights to the
   unweighted velocity integral used by the solver.
2. Use `lambda=1/T` consistently.  The previous Poiseuille initialization used
   `1/(2T)`.
3. Use case-dependent Prandtl numbers: 1 for Poiseuille and 2/3 for RT.
4. Apply `(ax,ay)=(phi,0)` only for Poiseuille and the radial acceleration only
   for RT; project that vector into the face normal/tangent frame.
5. Use the force-shifted interface velocity `vn-0.5*fn*dt`.
6. Correct the equilibrium time-slope moment call, the missing density factor
   in its tangential term, and all required quadrature weights.
7. Correct the reduced-B force term: the kinetic-energy factor multiplies
   `B0`, not `H0`.
8. Correct the sign of the horizontal-face tangential slope.
9. Apply specular and diffuse conditions only to the incoming velocity
   half-range and verify zero boundary mass flux.
10. Update every RT physical cell.  The legacy routine silently skipped the
    last row and column and then overwrote them.
11. Use `atan2(y,x)` for RT instead of a one-argument arctangent.
12. Use periodic ghost values for the Poiseuille streamwise reconstruction.
13. Normalize all steady residuals by the total-mass scale.  The legacy
    component-wise denominator becomes singular for the physically zero
    Poiseuille wall-normal momentum and can prevent an otherwise converged
    solution from ever satisfying the stopping condition.
14. Restore the case-dependent transport normalization.  Poiseuille uses
    `mu_ref=0.5*sqrt(pi)*Kn` and `tau=mu_ref/rho`, as in its original code;
    RT retains the temperature-dependent VHS collision time.  Reusing the RT
    formula in Poiseuille made its effective collision time about 10--20%
    too large and systematically underpredicted the paper's U and T profiles.

## Automated checks

`tests/run_tests.ps1` is compiled with
`-fcheck=all -ffpe-trap=invalid,zero,overflow`.

- The corrected two-dimensional symmetric-tensor implementation is tested
  independently with a manufactured rank-three Hermite expansion.  Its
  pointwise force is compared with the analytic expression using the direct
  \(d_{r,s}/(r!s!)\) contraction.
- The normalized \(b^{(n+1)}\) implementation forms only the independent
  components
  \[
  b^{(n+1)}_{r,s}=
  \frac{r a_xd^{(n)}_{r-1,s}+s a_yd^{(n)}_{r,s-1}}{n+1},
  \qquad r+s=n+1,
  \]
  and is compared pointwise with the direct-\(d\) result for the manufactured
  expansion and two shifted-Maxwellian acceleration directions.
- The new and legacy Hermite force arrays are compared for a shifted
  Maxwellian under Cartesian and oblique accelerations.
- GH16, GH18 and GH20: ordered/symmetric nodes and Maxwellian mass, first and
  second moments.
- Cartesian and radial Hermite force: mass, momentum and force-work moments.
- Left/bottom and right/top specular faces: half-range reflection and zero mass
  flux.
- Lower/upper diffuse Poiseuille walls: zero mass flux.
- Case-aware acceleration and periodic ghost-cell copies.

All regression checks pass.

## Numerical results

### Rayleigh--Taylor

Both grids reached the paper's `Kn=0.01` output time `t=0.14`.

| grid | steps | relative mass change | final rho range | final T range |
|---:|---:|---:|---:|---:|
| 60 x 60 | 236 | -1.36e-15 | 0.0023556 to 0.4849893 | 0.5425271 to 1.1196367 |
| 120 x 120 | 470 | 6.20e-15 | 0.0022792 to 0.4928833 | 0.5425047 to 1.1194282 |

After averaging each 2 x 2 block of the fine density onto the coarse grid, the
60/120 relative L2 differences are:

| time | relative L2 | absolute Linf |
|---:|---:|---:|
| 0.04 | 7.34e-3 | 1.09e-2 |
| 0.08 | 4.22e-3 | 4.00e-3 |
| 0.14 | 2.65e-3 | 1.82e-3 |

The four density snapshots reproduce the interface location and evolution of
the `Kn=0.01` row of the paper's RT figure.  This is a visual comparison because
the paper does not provide machine-readable reference fields.

### Poiseuille

#### Correct symmetric-tensor implementation

The legacy `force_type='hermite'` and the new
`force_type='hermite_symmetric'` were rebuilt from the same source and run in
separate directories with the same `3 x 100`, GH16, `Nm=7`, `t=0.2`
configuration.  Both runs completed 248 updates.  Their final `macro248.dat`
fields differ only at roundoff level:

| field | relative L2 difference | absolute Linf difference |
|---|---:|---:|
| rho | 3.42e-16 | 8.88e-16 |
| U | 1.34e-15 | 7.29e-17 |
| V | 9.55e-11 | 3.76e-16 |
| T | 6.04e-16 | 2.22e-15 |
| p | 5.91e-16 | 1.33e-15 |
| qx | 5.66e-13 | 1.14e-16 |
| qy | 6.52e-11 | 3.05e-16 |

The relatively larger normalized differences in `V` and `qy` are caused by
their near-zero reference norms; their absolute differences remain of order
`1e-16`.  The symmetric-tensor run has zero streamwise variation of `U`, a
mid-plane `U` symmetry error of `7.29e-17`, and maximum `|V|=3.54e-6`.
This verifies numerical equivalence of the two force operators for the tested
trajectory.  It does not replace the separate steady-state comparison with
the paper.

These short `t=0.2` trajectories are implementation checks only and must not
be compared with the paper's steady profiles.

#### Normalized symmetric \(b^{(n+1)}\) implementation

The direct-\(d\) case and `force_type='hermite_symmetric_b'` were rebuilt from
the same source and run with the same `3 x 100`, GH16, `Nm=7`, `t=0.2`
configuration.  Both completed 248 updates.  The symmetric-\(b\) run had a
relative mass change of `-4.16e-15`.  Its final field differences from the
direct-\(d\) run were:

| field | relative L2 difference | absolute Linf difference |
|---|---:|---:|
| rho | 3.27e-16 | 8.88e-16 |
| U | 1.25e-15 | 6.59e-17 |
| V | 8.28e-11 | 3.81e-16 |
| T | 6.41e-16 | 1.78e-15 |
| p | 5.46e-16 | 6.66e-16 |
| qx | 4.90e-13 | 1.09e-16 |
| qy | 7.71e-11 | 2.76e-16 |

The larger relative values for `V`, `qx`, and `qy` result from division by
near-zero reference norms.  Every absolute field difference is at roundoff
level.  Thus the explicit normalized symmetric-\(b^{(n+1)}\) contraction is
numerically equivalent to the direct-\(d^{(n)}\) contraction in this tested
Poiseuille trajectory.

The corrected 3 x 100 run stopped automatically when all four mass-scaled
residuals were below `1e-7`:

| quantity | result |
|---|---:|
| updates / final time | 48,315 / 61.4840 |
| relative mass change | -1.39e-14 |
| rho range | 0.952800 to 1.233240 |
| U range | 0.0888655 to 1.161727 |
| T range | 1.063461 to 1.373169 |
| maximum absolute V | 4.29e-4 |
| maximum x variation of U | 0 |
| U midplane symmetry error | 1.44e-15 |

The paper figure gives approximately `rho=0.95--1.22`,
`U=0.10--1.16`, and `T=1.07--1.375`.  The computed GH16, Nm=7 curves
have the same shape and agree with those plotted ranges to visual plotting
precision.  Exact curve norms cannot be given because the article contains no
machine-readable reference profile.

An intermediate diagnostic also isolated the transport-model error: at 10,000
steps the corrected Poiseuille model gave `Umax=0.949`, while reusing the RT
collision model gave `Umax=0.868`.  At 30,000 steps the incorrect model was
settling near `Umax=1.047` and `Tmax=1.298`, inconsistent with the paper; the
case-dependent model reached the reported paper ranges.

## Final runtime checks

The complete source was rebuilt after all edits.  All unit/regression checks
passed again.  Full main-program smoke cases were additionally compiled with
`-fcheck=all -ffpe-trap=invalid,zero,overflow`:

- Poiseuille 3 x 8: relative mass change 0, no runtime exception.
- RT 8 x 8: relative mass change `6.14e-16`, no runtime exception.

## Reproduction

```powershell
cd D:\ugks_hermite_force\ugks2d
.\build.ps1
.\tests\run_tests.ps1

mkdir results\poiseuille_steady
cd results\poiseuille_steady
..\..\ugks2d.exe ..\..\cases\poiseuille_steady.namelist

cd ..\..
python plot_hermite_validation.py
```

The plotting script is intentionally a plain Matplotlib script.  It writes
static PNG, PDF and EPS figures to `figures`.

## Scope and limitations

- Only Gauss--Hermite velocity quadrature and the Hermite force expansion were
  exercised.
- Passing these cases does not validate all Knudsen numbers or all boundary
  configurations.
- The paper lists a non-square RT mesh in one table, while its supplied RT
  source constructs a 120 x 120-cell square grid.  The square source geometry,
  the plotted quarter-circle domain, and the user's requested 60 x 60 then
  120 x 120 sequence were used here.
- Agreement with a plotted paper curve/contour is not an independent physical
  validation; raw reference arrays would be required for a strict pointwise
  error norm against the article.
