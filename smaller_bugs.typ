= Smaller Radiation Bug Fixes

#set page(margin: 2cm)
#set text(size: 11pt)

== Smaller Bug 1: Axis boundary uses wrong `isTwoPi` flag

*File:* `src/radiation/FLD.cpp`, `FluxLimitedDiffusion::EnforceBoundary`

`FluxLimitedDiffusion::Init` correctly determines whether the domain spans the full $2\pi$ in $\phi$ and stores the result in `rad_boundary.isTwoPi`. However `EnforceBoundary` read the class member `this->isTwoPi`, which is a separate `bool` that is never set, so it is always `false`.

In a full-$2\pi$ domain the axis ghost cells across the pole should be filled by mirroring at $\phi + \pi$. With `isTwoPi = false` they were mirrored at the same $\phi$, corrupting the flux-limiter gradient and the radiation pressure in the cells adjacent to the pole.

*Fix:* Pass `rad_boundary.isTwoPi` instead of `this->isTwoPi` to `ApplyAxisBoundary`.

#pagebreak()
