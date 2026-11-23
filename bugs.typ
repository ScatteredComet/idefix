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
