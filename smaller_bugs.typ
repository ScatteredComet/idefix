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
