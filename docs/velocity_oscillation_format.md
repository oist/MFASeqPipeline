# Velocity oscillation file format

Produced by the `velocity_oscillation` step of [mfa_pipeline](mfa_pipeline.md)
(`steps/fit_velocity_oscillation.R`). The step fits a model of the replication
fork speed to the abundance profile, assuming a circular genome with a single
origin.

Tab-separated text with a header line and four rows, one per model regime:

| `regime` | Model |
|---|---|
| `vc.Dz` | constant speed, no diffusion (`delta = 0`, `D = 0`) |
| `vosc.Dz` | oscillating speed, no diffusion (`D = 0`) |
| `vc.Dnz` | constant speed with diffusion (`delta = 0`) |
| `vosc.Dnz` | oscillating speed with diffusion |

## Columns

| Column | Description |
|---|---|
| `model` | Model family, the `velocity_oscillation_model` setting: `vt` (speed oscillates in time), `vx` (speed oscillates with position) or `vd` (speed oscillates with distance travelled) |
| `regime` | See above |
| `n` | Number of windows used: those with a finite `coverage_ratio` and `ratio_mappable = 1` |
| `df` | Number of free parameters in the regime |
| `optimal` | `TRUE` for the regime with the lowest AIC (ties broken by fewer `df`), otherwise `FALSE` |
| `L` | Genome length in bp |
| `O` | Origin position in bp |
| `lambda` | Growth rate in 1/s, i.e. `growthrate_perhour / 3600` |
| `D` | Diffusion coefficient of the fork position |
| `v0` | Mean fork speed, in bp/s |
| `delta` | Relative amplitude of the speed oscillation (0 ≤ `delta` ≤ 1) |
| `omega` | Angular frequency of the oscillation: rad/s for `vt`, rad/bp for `vx` and `vd` |
| `phi` | Phase of the oscillation, in rad |
| `residual` | Weighted residual sum of squares |
| `logl` | Log-likelihood, assuming normally distributed errors with the standard errors from the abundance file |
| `aic` | Akaike information criterion, `2 df - 2 logl` |

In the `vt` model, the fork speed at time `t` after initiation is
`v(t) = v0 * (1 + delta * cos(omega * t + phi))`. The `vx` and `vd` models use
the same formula with distance from the origin, or distance travelled, in place
of `t`.
