# In Vitro Model Fits

Posterior and prior summaries for MCMC fits to Penn *et al.* in vitro A549 data (PB1 full-genome, DVG, IFN-β mRNA).

## Files

- **posterior_table.md** — Posterior medians and 90% credible intervals [5%, 95%] for base model and our model across low/med/high DVG groups
- **prior_table.md** — Prior specifications and inference settings

## Data source

Penn *et al.* (2022) PB1 full-genome area, PB1 DVG area, and IFN-β mRNA from A549 cells infected with 6:2 Tky/05 (low DVG), 7:1 Tky/05 LOW (med DVG), or 7:1 Tky/05 HIGH (high DVG) at 0, 2, 6, 24 h.

## Inference

- **Sampler:** emcee (Python)
- **Walkers:** 32
- **Burn-in:** 3000 steps
- **Production:** 20,000 steps per walker

Raw MCMC samples and fitting scripts: `0204_in_vitro/` in the re-evalu_IAV_DIP_kinetics repository.
