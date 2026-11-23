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
