# Julia code: Fisher Information and MBAM analysis

This folder contains the Julia code for the Fisher Information Matrix (FIM)
analysis and the Manifold Boundary Approximation Method (MBAM) reduction of the
infinite-bus synchronous generator in

> N. Evangelou, A. M. Stankovic, I. G. Kevrekidis, and M. K. Transtrum,
> "Model-based and data-driven approaches to parameter identifiability:
> A power systems case study."

Author: Mark K. Transtrum.

## Layout

```
Project.toml, Manifest.toml   Julia environment (Julia 1.10.4)
src/
  mbam_longtimes.jl           Driver: FIM and the five MBAM iterations (11 -> 6 parameters)
  pvalues.jl                  Nominal parameters and constants
  SG11.jl ... SG06.jl         Parameterizations for each step of the reduction
  ic.jl, ic08.jl, ic07.jl     Initial conditions (and consistent derivatives for the DAE solver)
  rhs.jl                      Model equations (ODE right-hand side and DAE residual)
  obs.jl                      Observation function (full state)
data/                         Outputs of mbam_longtimes.jl used for the figures
figures/                      Scripts that produce Figs. 2-7 of the paper from data/
check/fim_check.py            Independent FIM check in Python (NumPy/SciPy only); see "Known issue"
```

## Model and parameterization

`rhs.jl` implements the full model of Sec. III of the paper. Parameters are the
"independent" parameters of Table IV (`deltax1`-`deltax5`, `deltaTd`, `deltaTq`,
`x2d`, `T2d0`, `T2q0`, `H`, `D`). The code uses `1`/`2` for single/double
primes, e.g. `x2d` is $x''_d$ and `T2d0` is $T''_{d0}$. All varied parameters
are log-transformed (`exp` in the `@parameterspace` definitions), so
sensitivities are with respect to log-parameters. `deltax5 = 0` is held fixed,
leaving 11 parameters.

The DAE residual `res` in `rhs.jl` rescales the time constants with
`xform(T) = T/(1+T)` (or `1/(1+1/T)` for `T > 1`), so that the singular limits
`T -> 0` can be taken without dividing by zero.

Each `SGnn.jl` file defines the model after a step of the reduction:

| File | Parameters | Limit applied to reach the next model |
|------|-----------|-------|
| `SG11.jl`  | 11 | $D \to 0$ |
| `SG10b.jl` | 10 | $H \to 0$ |
| `SG09.jl`  |  9 | $T''_{d0} \to 0$ |
| `SG08.jl`  |  8 | $T''_{q0} \to 0$ |
| `SG07.jl`  |  7 | $\delta x_1 \to 0$ |
| `SG06.jl`  |  6 | (reduced model) |

The limit $H \to 0$ is approximated numerically by fixing `H = 1e-1` in
`SG09.jl` through `SG06.jl`, rather than setting it to zero. This is a
deliberate numerical stand-in that avoids reformulating the swing equation as
an algebraic constraint. The limit $D \to 0$ is taken exactly.

### Values used

- Parameters: Table III of the paper.
- Constants (`pvalues.jl`): `freq = 60`, `omega0 = 1` (per unit),
  `omegab = 2*pi*freq`, `Pm = 0.7`, `V = 1.09`, and `vf0 = 4.14676`,
  `theta = -0.0382548` (taken from a PSAT simulation of the system).
- Initial condition (`ic.jl`): `delta = 0.5`, `omega = 0.98`, `e'_q = 2.13`,
  `e'_d = 0.02`, `e''_q = 1.98`, `e''_d = 0.02`.
- These differ slightly from Tables I and II of the paper, which list
  `vf0 = 4.2`, `theta = 0` and `e''_q = 1.93`. The difference in `e''_q(0)`
  has a negligible effect on the observed window `t = 3:5`. The difference in
  `vf0` and `theta` shifts the operating point (by up to 0.03 in `delta`) but
  leaves the FIM spectrum and the identifiable parameters essentially unchanged.
- Observations: all six states at `t = 3:0.02:5`, giving a 606 x 11 Jacobian.
- Solver tolerances: `abstol = reltol = 1e-12`.

## Known issue: $x''_d$ is held fixed in the model equations

The parameter $\delta x_5 = x''_q - x''_d$ is zero at the nominal values, so it
lies on a boundary of parameter space. It is therefore fixed at zero, leaving 11
parameters, and varying `x2d` is meant to vary $x''_q = x''_d$ together.

In this code, however, the function `ps()` in each `SGnn.jl` file returns the
physical parameters *without* `x2d`. The model functions in `rhs.jl` and the
`ic*.jl` files therefore use the global value `x2d = 0.48` from `pvalues.jl`.
As a result, varying the `x2d` parameter changes $x''_q$ (and, through the
chain of differences, $x'_d$, $x'_q$, $x_q$, $x_d$), but not the $x''_d$ that
appears in $i_d$ and in the $e''_q$ equation. In effect, the parameter labeled
$x''_d$ in this analysis acts as $x''_q$, with $x''_d$ fixed.

Nominal trajectories (Figs. 2 and 6) are unaffected, because the two values
agree at the nominal point. The sensitivities, and therefore the FIM spectra
and geodesics, are affected. The code is left as it was run for the paper,
because the packages it depends on are no longer available to rerun it.

`check/fim_check.py` reimplements the model with only NumPy and SciPy. It
reproduces the saved Jacobian `data/SG11_jacobian.txt` to a relative difference
of $6\times10^{-9}$ when it replicates the behavior above, then recomputes the
FIM with $x''_d$ varied consistently:

| | FIM eigenvalues |
|---|---|
| As used in the paper | 140, 13.2, 5.37, 0.289, 0.209, 8.2e-3, 1.3e-4, 2.9e-5, 6.5e-6, 8.2e-8, 1.2e-9 |
| Corrected | 145, 10.0, 5.49, 3.75, 0.221, 9.6e-3, 2.9e-4, 3.0e-5, 6.5e-6, 8.2e-8, 1.3e-9 |

The conclusions are unchanged. In both cases six eigenvalues are of order
$10^{-2}$ or larger, followed by a gap, and the six stiffest directions are
spanned by the same six parameters ($\delta x_2$, $\delta x_3$, $\delta x_4$,
$x''_d$, $\delta T_d$, $\delta T_q$, each with weight at least 0.99). The
smallest eigenvalue, which sets the length of the first geodesic ($D \to 0$),
is also essentially unchanged.

## Running

From `src/`:

```julia
include("mbam_longtimes.jl")
```

This regenerates the files in `data/` and writes diagnostic plots to
`src/MBAM/`. The figure scripts in `figures/` read from `data/` and require
`CairoMakie` and `LaTeXStrings`. As written, they activate a shared environment
named `Plots` (`Pkg.activate("Plots"; shared = true)`); replace that line with
any environment that provides those packages.

### Dependencies

`Models`, `ParametricModels`, and `Geometry` are internal research packages
from the Transtrum group. They are **not** publicly available and are not in
Julia's General registry (the `Models` package in the General registry is an
unrelated package), so this code cannot be run as-is outside the group.
`Manifest.toml` records the exact versions used. The code is provided to
document the analysis: the model equations, parameterizations, initial
conditions, and the sequence of MBAM limits. The saved outputs in `data/` are
sufficient to reproduce the figures with only public packages.

The packages are used for: the parameter-space definitions
(`@parameterspace`, `parameter_transforms`, `xvalues`), solving the ODE/DAE
models (`solve_ode`, `solve_dae`), forward-mode sensitivities
(`model.jacobian`), second directional derivatives (`model.Avv`), and
geodesic integration (`Geometry.geodesic`), following
Transtrum, Machta and Sethna, *Phys. Rev. E* **83**, 036701 (2011).
