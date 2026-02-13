# Supplementary Table S3: ODE Model Comparison (V+F+D models, n=60 per group)

BIC for low-D₀, intermediate-D₀, and high-D₀ inoculum groups. **Our model achieves the lowest BIC in all three groups.**

**BIC calculation:** BIC = −2·elpd + ln(n)·k, where elpd is the sum over observations of the posterior mean log-likelihood (computed from Stan MCMC draws), n = 60 for V+F+D models, and k is the number of free parameters (excluding initial conditions V₀, F₀, D₀ and nuisance parameters). Lower BIC indicates better fit with parsimony.

Full Stan fits are available at: https://github.com/yimei-li/within-host-dip-iav-ode

---

## Our model (best)

| # | Full ODE system (T, I, V, D, F) | Comments / Stan inference | BIC low | BIC med | BIC high |
|---|--------------------------------|--------------------------|---------|---------|---------|
| 0 | **T:** dT/dt = −β T V<br>**I:** dI/dt = β T V − δ I<br>**V:** dV/dt = (1−d) p_eff I − c V, with p_eff = p̂/(1+ε₂F)<br>**D:** dD/dt = d p_eff I − h D<br>**F:** dF/dt = s(D) I − α F, with s(D) = s₀ + s_d·D/(1+θ_s·D) | Stable convergence. No divergent transitions. | **254.2** | **249.6** | **235.1** |

---

## Alternatives

| # | Full ODE system (T, I, V, D, F) | Comments / Stan inference | BIC low | BIC med | BIC high |
|---|--------------------------------|--------------------------|---------|---------|---------|
| 1 | **T:** dT/dt = −k_eff T V, with k_eff = k̂/(1+ε₃D)<br>**I:** dI/dt = k_eff T V − δ I<br>**V:** dV/dt = (1−d) p_eff I − c V, with p_eff = p̂/(1+ε₂F)<br>**D:** dD/dt = d p_eff I − h D<br>**F:** dF/dt = s(D) I − α F, with s(D) = s₀ + s_d·D/(1+θ_s·D) | ε₃ narrow (ca. 10⁻³); weak DIP-inhibition. | 272.3 | 266.6 | 245.2 |
| 2 | **T:** dT/dt = −β_eff T V, with β_eff = β̂/(1+ε₃D)<br>**I₁:** dI₁/dt = β_eff T V − k_e I₁<br>**I₂:** dI₂/dt = k_e I₁ − δ I₂<br>**V:** dV/dt = (1−d) p I₂ − c V<br>**D:** dD/dt = d p I₂ − h D<br>**F:** dF/dt = s(D) I₂ − α F, with s(D) = s₀ + s_d·D/(1+θ_s·D) | Same ε₃ pattern. Eclipse phase I₁→I₂. | 273.2 | 267.3 | 245.6 |
| 3 | **T:** dT/dt = −k_eff T V, with k_eff = k̂/(1+ε₁F)<br>**I:** dI/dt = k_eff T V − δ I<br>**V:** dV/dt = (1−d) p_eff I − c V, with p_eff = p̂/(1+ε₂F)<br>**D:** dD/dt = d p_eff I − h D<br>**F:** dF/dt = s(D) I − α F, with s(D) = s₀ + s_d·D/(1+θ_s·D) | k̂, ε₁ correlated. | 271.0 | 266.2 | 253.1 |
| 4 | **T:** dT/dt = −k_eff T V, with k_eff = k̂/(1+ε₁F+ε₃D)<br>**I:** dI/dt = k_eff T V − δ I<br>**V:** dV/dt = (1−d) p_eff I − c V, with p_eff = p̂/(1+ε₂F)<br>**D:** dD/dt = d p_eff I − h D<br>**F:** dF/dt = s(D) I − α F, with s(D) = s₀ + s_d·D/(1+θ_s·D) | Overparameterised. | 280.2 | 275.0 | 253.7 |
| 5 | **T:** dT/dt = −β T V<br>**I₁:** dI₁/dt = β T V − k_e I₁<br>**I₂:** dI₂/dt = k_e I₁ − δ I₂<br>**V:** dV/dt = (1−d) p I₂ − c V<br>**D:** dD/dt = d p I₂ − h D<br>**F:** dF/dt = s(D) I₂ − α F, with s(D) = s₀ + s_d·D/(1+θ_s·D) | k_e identifiable; no gain. Eclipse I₁→I₂; no ε. | 264.8 | 267.3 | 254.9 |
| 6 | **T:** dT/dt = −β_eff T V, with β_eff = β̂/(1+ε₃D)<br>**I₁:** dI₁/dt = β_eff T V − k_eff I₁, with k_eff = k_e/(1+ε₁F)<br>**I₂:** dI₂/dt = k_eff I₁ − δ I₂<br>**V:** dV/dt = (1−d) p_eff I₂ − c V, with p_eff = p/(1+ε₂F)<br>**D:** dD/dt = d p_eff I₂ − h D<br>**F:** dF/dt = s(D) I₂ − α F, with s(D) = s₀ + s_d·D/(1+θ_s·D) | Overparameterised. β_eff, k_eff(F), p_eff(F); eclipse. | 288.8 | 283.1 | 262.0 |

---

## Model key

- **T:** target cells  
- **I, I₁, I₂:** infected cells (single compartment or eclipse phase I₁→I₂)  
- **V:** infectious virus  
- **D:** DIP/DVG signal  
- **F:** IFN-α  

Parameters: β (infection rate), δ (infected-cell loss), p̂ (virus production), c (viral clearance), h (DIP clearance), α (IFN clearance), d (DIP fraction), s₀, s_d, θ_s (IFN induction), ε₁, ε₂, ε₃ (inhibition strengths), k_e (eclipse rate).
