# In Vitro Prior Specifications

Uniform priors on the specified ranges. Inference: emcee MCMC (32 walkers, 3000 burn-in, 20,000 production steps per walker).

---

## Base model (TIVF)

| Parameter | Description | Prior range [min, max] |
|-----------|-------------|------------------------|
| β | Infection rate | [1e-12, 0.001] |
| δ | Infected-cell death rate (day⁻¹) | [0.01, 10] |
| p̂ | Virus production rate | [0.1, 100] |
| ε₂ | IFN inhibition of virus production | [0.01, 0.95] |
| c | Virus clearance rate (day⁻¹) | [0.01, 3] |
| s | IFN production rate (base model) | [1e-8, 1] |
| α | IFN clearance rate (day⁻¹) | [0.01, 30] |
| V₀ | Initial virus | [0.01, 6624] (low), [0.01, 4533] (med), [0.01, 459] (high) |
| F₀ | Initial IFN | [0.0001, 5.43] (low), [0.0001, 26.1] (med), [0.0001, 130] (high) |
| σ_V | Virus observation error (log scale) | [0.1, 2] |
| σ_F | IFN observation error (log scale) | [0.1, 2] |

V₀ and F₀ upper bounds are data-dependent (max observed × 1.2 for each group).

---

## Our model (TIVFD)

| Parameter | Description | Prior range [min, max] |
|-----------|-------------|------------------------|
| β | Infection rate | [1e-12, 0.001] |
| δ | Infected-cell death rate (day⁻¹) | [0.01, 10] |
| p̂ | Virus production rate | [0.1, 100] |
| ε₂ | IFN inhibition of virus production | [0.01, 0.95] |
| c | Virus clearance rate (day⁻¹) | [0.01, 3] |
| d | DIP production fraction | [1e-6, 0.999] |
| h | DIP clearance rate (day⁻¹) | [0.01, 3] |
| s_d | DIP-driven IFN induction | [0.0001, 20] |
| θ_s | DIP–IFN saturation scale | [1, 100000] |
| α | IFN clearance rate (day⁻¹) | [0.01, 30] |
| V₀ | Initial virus | [0.01, 6624] (low), [0.01, 4533] (med), [0.01, 459] (high) |
| D₀ | Initial DIP | [0.01, 3500] (low), [0.01, 6225] (med), [0.01, 10983] (high) |
| F₀ | Initial IFN | [0.0001, 5.43] (low), [0.0001, 26.1] (med), [0.0001, 130] (high) |
| σ_V | Virus observation error (log scale) | [0.1, 2] |
| σ_D | DIP observation error (log scale) | [0.1, 2] |
| σ_F | IFN observation error (log scale) | [0.1, 2] |

V₀, D₀, and F₀ upper bounds are data-dependent per DVG group.

---

## Inference settings

- **Sampler:** emcee `EnsembleSampler`
- **Walkers:** 32
- **Burn-in:** 3000 steps
- **Production:** 20,000 steps per walker
- **Total samples:** 640,000 (32 × 20,000)

Full scripts and raw MCMC samples: see `0204_in_vitro/` in the re-evalu_IAV_DIP_kinetics repository.
