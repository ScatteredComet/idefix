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
