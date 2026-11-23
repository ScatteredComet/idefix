= Radiation Bug Fixes

#set page(margin: 2cm)
#set text(size: 11pt)

== Bug 1: Long-characteristic irradiation counts each cell's own optical depth twice

*File:* `src/radiation/irradiation.cpp`, `ComputeHeating`

*Problem:* `Column::ComputeColumn` returns an _inclusive_ prefix sum, so `column(s,k,j,i)` already contains the contribution $\rho_s \, \Delta r$ of the current cell. The heating formula computes

$
H = F_0 \, e^{-\tau} \, (1 - e^{-\Delta\tau})
$

where $\tau$ is the optical depth upstream of the cell and $\Delta\tau = \kappa_s \rho_s \Delta r$ is the cell's own optical depth. With the inclusive column, $\tau$ was set equal to the cumulative sum through and including the cell, so the cell's own $\Delta\tau$ appeared both in the $e^{-\tau}$ factor and in $(1-e^{-\Delta\tau})$. In the optically thick limit this suppresses absorption by a factor $e^{-\Delta\tau}$.

*Fix:* Use the exclusive upstream optical depth $\tau = (\text{column} - \Delta\text{col}) \cdot \kappa$.

```cpp
// Before:
tau += column(s, k, j, i) * kappa(s, l);
// After:
tau += (column(s, k, j, i) - dcol(s, k, j, i)) * kappa(s, l);
```

#pagebreak()

== Bug 2: Two-temperature FLD coupling coefficient

*File:* `src/radiation/FLD.cpp`, `TwoTemperatureFLD::FillMatrixCouplingTerms`

The FLD equation for the radiation energy density $E$ in the two-temperature scheme is solved implicitly. After eliminating the gas temperature, the diagonal coefficient of the matrix for $E$ is

$
d_{ii} = 1 + \frac{y}{1 + 4x}
$

where $y = \kappa_P \rho \, c \, \Delta t$ and $x = a_R T^3 y / (\rho c_V)$. The code had

```cpp
// Before (wrong):
diag(k, j, i) = (1 + radY(0, k, j, i)) / (1 + 4 * radX(0, k, j, i));
// After (correct):
diag(k, j, i) = 1 + radY(0, k, j, i) / (1 + 4 * radX(0, k, j, i));
```

The incorrect formula divided both $1$ and $y$ by $(1+4x)$, shifting the diagonal below 1 in the optically thick limit and making the implicit system under-constrained.

#pagebreak()

== Bug 3: Two-temperature FLD irradiation of dust missing density weighting

*File:* `src/radiation/FLD.cpp`, `TwoTemperatureFLD::UpdatePressure`

The irradiation heating rate $S_s$ (per unit mass of species $s$) contributes to the gas temperature update as

$
\delta T_{\text{gas}} = \frac{\Delta t}{c_V} \sum_s S_s \frac{\rho_s}{\rho_{\text{gas}}}
$

because the total energy deposited per unit volume in species $s$ is $S_s \rho_s$, which must be divided by $\rho_{\text{gas}} c_V$ to get a gas temperature increment. The code had

```cpp
// Before (wrong):
for (int s = 0; s < num_species; s++) num += Sirrad(s, k, j, i) * dt / cV;
```

which omitted the $\rho_s/\rho_{\text{gas}}$ factor for dust species ($s > 0$), overheating the gas by $1/\epsilon$ where $\epsilon$ is the dust-to-gas ratio. The corresponding `FillUtils` path already contained the correct weighting.

*Fix:* Pre-compute the density-weighted sum before the main kernel and use it in the numerator.

#pagebreak()

== Bug 4: FLD radiation-pressure `Vc` alias pointed to `Uc`

*File:* `src/radiation/FLD.cpp`, `FluxLimitedDiffusion::_ComputeRadiationPressureSourceTerm`

In the branch that updates primitive variables (`update_Vc = true`), the code had

```cpp
// Before:
IdefixArray4D<real> Vc = data->hydro->Uc;  // wrong alias
```

so `Vc` and `Uc` both pointed to the conservative array. The momentum kick was therefore written into `Uc`, which is then immediately overwritten by the subsequent `PrimToCons` call in the time integrator, silently discarding the radiation-pressure kick on the gas.

*Fix:* `IdefixArray4D<real> Vc = data->hydro->Vc;`

#pagebreak()

== Bug 5: Radiation-pressure momentum kick missing energy update

*Files:* `src/radiation/FLD.cpp`, `src/radiation/irradiation.cpp`, `src/radiation/short_char.cpp`

All three `_ComputeRadiationPressureSourceTerm` implementations (FLD, long-char, short-char) updated the conservative momentum density $\mathbf{m} = \rho\mathbf{v}$, but never updated the total energy $E = \frac{1}{2}\rho v^2 + e$. The work done by the radiation force on the fluid is

$
\frac{\partial E}{\partial t} = \mathbf{v} \cdot \mathbf{f}_{\rm rad}
$

Omitting this means the kinetic energy increase from the momentum kick comes at the expense of the internal energy, which is physically wrong.

*Fix:* After adding momentum kick $\delta\mathbf{m}$, add the work term

```cpp
Uc(ENG, k, j, i) += Vc(VX1, k, j, i) * dp1
                  + Vc(VX2, k, j, i) * dp2  // if DIMENSIONS >= 2
                  + Vc(VX3, k, j, i) * dp3; // if DIMENSIONS == 3
```

guarded by `#if HAVE_ENERGY == 1`.

#pagebreak()

== Bug 6: Dusty-wind test applies radiation pressure twice

*File:* `test/Radiation/dynamic/dusty-wind/setup.cpp`

`RadiationPressureSourceTerm` is enrolled as a user source term via `EnrollUserSourceTerm`, so the time integrator already calls it every RK stage. The same function was then called again inside `UpdateSurface` (registered as `UserStepLast`), doubling the radiation-pressure kick for all species every time-step.

Note: bug 4 happened to conceal the double kick on the gas (the second `update_Vc = true` kick landed in `Uc` and was then discarded by `PrimToCons`), but dust received an unambiguous double kick.

*Fix:* Remove the redundant `ComputeRadiationPressureSourceTerm` block from `UpdateSurface`.

#pagebreak()

== Bug 7: Levermore–Pomraning limiter diverges at small $R$

*File:* `src/radiation/FLD.cpp`, `RadiationFluxLimiterLevermorePomraning1981`

The limiter is

$
\lambda(R) = \frac{1}{R}\left(\coth R - \frac{1}{R}\right)
$

whose correct limit as $R \to 0$ is $\lambda \to 1/3$. The code guarded against division by zero with

```cpp
real one_over_R = 1.0 / (R + 1e-15);
return one_over_R * (1.0 / tanh(R) - one_over_R);
```

For $R \lesssim 10^{-6}$, `1/tanh(R)` and `1/(R+1e-15)` are both $\approx R^{-1}$ and their difference underflows to zero through catastrophic cancellation, leaving $\lambda \approx 0$. For $R \lesssim 10^{-15}$, `1/(R+1e-15)` ≈ $10^{15}$, causing $\lambda \gg 1$. Optically thick interiors routinely fall in this range.

*Fix:* Use the Taylor series $\lambda \approx 1/3 - R^2/45$ for $R < 0.01$.

#pagebreak()

== Bug 8: Long-characteristic boundary column has wrong units

*File:* `src/radiation/irradiation.cpp`, `LongCharIrradiation::ComputeIrradiation`

`boundaryColumn(s, k, j)` stores the column density (mass per unit area) from outside the domain to the inner radial boundary. It should be added as a constant offset to the cumulative column integral. Instead it was added to `dcol` (the volumetric mass density) before the scan:

```cpp
// Before:
dcol(s,k,j,i) = Vc(RHO, k, j, i);
if (rank == 0 && i == ibeg) dcol(s,k,j,i) += boundaryColumn(s, k, j);
column_sum->ComputeColumn(dcol, s);   // integrates dcol * dr
column(s,k,j,i) = col_result(k,j,i);
```

`ComputeColumn` integrates $\rho\, \mathrm{d}r$, so `boundaryColumn` was multiplied by the first cell's path length $\Delta r_0$ instead of being taken as-is. For a cell with $\Delta r_0 \neq 1$, every downstream cell received an incorrect column offset.

*Fix:* Remove `boundaryColumn` from `dcol`; add it directly to `col_result` after the scan:

```cpp
column(s,k,j,i) = col_result(k,j,i) + boundaryColumn(s, k, j);
```

#pagebreak()

== Bug 9: Null `precond` view passed to user opacity with non-preconditioned solver

*File:* `src/radiation/FLD.cpp`, `FluxLimitedDiffusion::_ComputeRadiationPressureSourceTerm`

`StencilMatrix::precond` is only allocated when `enablePreconditioner()` is called, which only happens for PBICGSTAB, PCG, and PMINRES solvers. For the default BICGSTAB solver, `precond` is an empty (null) `IdefixArray3D`. The radiation-pressure routine passed it as a scratch array to the user opacity callback:

```cpp
// Before:
data->radiation->UserOpacityFuncCGS(data, species, FLD_matrix->precond, kappaR);
```

Any write to `FLD_matrix->precond` from within the opacity function would dereference a null pointer, causing a crash or silent memory corruption.

*Fix:* Allocate a local temporary array and pass that instead:

```cpp
IdefixArray3D<real> kappaP_scratch = IdefixArray3D<real>(
    "kappaP_scratch", np_tot[KDIR], np_tot[JDIR], np_tot[IDIR]);
data->radiation->UserOpacityFuncCGS(data, species, kappaP_scratch, kappaR);
```
