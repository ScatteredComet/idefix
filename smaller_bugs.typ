= Smaller Radiation Bug Fixes

#set page(margin: 2cm)
#set text(size: 11pt)

== Smaller Bug 1: Axis boundary uses wrong `isTwoPi` flag

*File:* `src/radiation/FLD.cpp`, `FluxLimitedDiffusion::EnforceBoundary`

`FluxLimitedDiffusion::Init` correctly determines whether the domain spans the full $2\pi$ in $\phi$ and stores the result in `rad_boundary.isTwoPi`. However `EnforceBoundary` read the class member `this->isTwoPi`, which is a separate `bool` that is never set, so it is always `false`.

In a full-$2\pi$ domain the axis ghost cells across the pole should be filled by mirroring at $\phi + \pi$. With `isTwoPi = false` they were mirrored at the same $\phi$, corrupting the flux-limiter gradient and the radiation pressure in the cells adjacent to the pole.

*Fix:* Pass `rad_boundary.isTwoPi` instead of `this->isTwoPi` to `ApplyAxisBoundary`.

#pagebreak()

== Smaller Bug 2: Stale hydro ghost cells before radiation solve

*File:* `src/timeIntegrator.cpp`, `TimeIntegrator::Cycle`

The stage loop ends with `ConsToPrim`, leaving interior ghost cells at the values set at the _beginning_ of the final stage. No `SetBoundaries` is called before `ComputeIrradiation` or `SolveSystem`. As a result:

- The diffusion coefficient $\nu = \lambda c / (\kappa_R \rho)$ is evaluated with ghost-cell densities and temperatures that are one stage behind. Two MPI ranks sharing a face therefore compute slightly different face coefficients, so the diffusion matrix is not exactly symmetric. This breaks CG and MINRES solvers (which rely on exact symmetry) and causes a small $O(\Delta t)$ non-conservation of radiation energy across MPI boundaries.

*Fix:* Call `data.SetBoundaries()` once before the radiation solve when radiation is active.

#pagebreak()

== Smaller Bug 3: Non-zero Neumann ghost-cell fill inconsistent with the matrix

*File:* `src/radiation/FLD.cpp`, `FluxLimitedDiffusion::EnforceBoundary`

The matrix enforces a Neumann boundary by folding the off-diagonal coupling into the diagonal and adding a correction

$
\text{rhs}_{i} \mathrel{+}= F \cdot M_{i,{i-1}} \cdot \Delta x_{\nu}
$

where $F$ is the specified outward flux and $\Delta x_\nu = \tfrac{1}{2}\bigl(dx_{i-1}/\nu_{i-1} + dx_i/\nu_i\bigr)$. This is equivalent to imposing $E_{\rm ghost} = E_{\rm ref} - F \cdot \Delta x_\nu$. The ghost-cell fill in `EnforceBoundary` ignored $F$ entirely and simply copied $E_{\rm ref}$, giving the wrong ghost-cell gradient for the flux limiter when $F \neq 0$.

All tests use $F = 0$, where both formulas agree.

*Fix:*
```cpp
localVar(k, j, i) = localVar(kref, jref, iref)
                 - val * dx_arr(dir_ref) / nu_arr(kref, jref, iref);
```
(approximating $dx_{\rm ghost} \approx dx_{\rm ref}$ and $\nu_{\rm ghost} \approx \nu_{\rm ref}$).

#pagebreak()

== Smaller Bug 4: `CYLINDRICAL` geometry silently uses a wrong metric

*File:* `src/radiation/FLD.cpp`, `_fill_FLD_matrix` and `_ComputeRadiationPressureSourceTerm`

Both routines compute the scale factor $h_2$ for the JDIR flux with

```cpp
#elif GEOMETRY == POLAR
  h2 = x1(i);   // correct: h_phi = R
#else
  h2 = x1(i);   // wrong for CYLINDRICAL: h_z = 1, not R
#endif
```

In `CYLINDRICAL` coordinates $(R,\,z)$ the $z$-direction has $h_2 = 1$ (there is no arc-length scale factor in $z$). The fallthrough `#else` branch assigns $h_2 = R$, introducing a spurious $1/R$ weighting in the $z$-diffusion flux and in the flux-limiter gradient. The same incorrect metric was used in the radiation-pressure term. There is no compile-time guard against this geometry.

*Fix:* Reject `CYLINDRICAL` geometry at startup with `IDEFIX_ERROR`.

#pagebreak()

== Smaller Bug 5: Dust irradiation without energy equation is silently discarded

*File:* `src/radiation/irradiation.cpp`, `Irradiation::UpdatePressure`

After computing the per-species irradiation heating array $S_s$ (energy per unit time per unit mass of species $s$), `UpdatePressure` applies $S_0$ to the gas pressure and applies $S_d$ ($d > 0$) to the dust tracer temperature only when `have_energy = true`. When dust has no energy equation, the absorbed starlight $S_d \rho_d \Delta t$ was simply dropped — neither stored in the dust nor transferred to the gas.

Physically, the absorbed energy must go somewhere. Without a dust energy variable, the nearest correct choice is to give it to the gas:

$
\delta P_{\rm gas} = S_d \,\rho_d \,\Delta t \,(\gamma - 1)
$

*Fix:* Add an `else` branch that routes the unabsorbed dust irradiation energy into the gas pressure.

#pagebreak()
