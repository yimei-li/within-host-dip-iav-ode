# Prior bounds and prior SD for Bayesian inference

Ranges were informed by MLE pilot fits and literature (Miao et al., Baccam et al., Pawelek et al.). Prior SD is the standard deviation on the log scale for log-normal priors (or on the logit scale for d); for σ_V, σ_D, σ_F we used exponential(0.5) priors.

## Shared by base and DIP-extended models

| Parameter | Bounds | Prior SD (log scale) |
|-----------|--------|----------------------|
| β | 10⁻⁸ – 10⁻³ | 3 |
| δ (1/day) | 0.2 – 8 | 0.7 |
| p̂ (virions/cell/day) | 0.1 – 1000 | 3.5 |
| ε / ε₂ | 10⁻⁴ – 10 | 3.5 |
| c (1/day) | 2 – 12 | 0.7 |
| α (1/day) | 1 – 15 | 0.7 |
| V₀ | 10⁻³ – 10⁶ | 2 |
| σ_V (obs.) | 0.1 – 5 | exp(0.5) |
| σ_F (obs.) | 0.1 – 1.5 | exp(0.5) |

## Base model only

| Parameter | Bounds | Prior SD (log scale) |
|-----------|--------|----------------------|
| s₀ | 10⁻⁴ – 1 | 2 |

## DIP-extended model only

| Parameter | Bounds | Prior SD (log scale) |
|-----------|--------|----------------------|
| d (DIP fraction) | (0, 1) | 2 (logit) |
| h (1/day) | 2 – 15 | 1.2 |
| s_d | 10⁻⁴ – 20 | 2 |
| θ_s | 1 – 10⁵ | 2 |
| D₀ | 1 – 10⁸ | 2 |
| σ_D | 0.1 – 5 | exp(0.5) |
