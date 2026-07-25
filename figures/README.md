# Figures

Pairwise posterior distributions for the with-DIP in vivo fit, one per DVG inoculum
group. 14 fitted parameters, log10 scale. Diagonal: 1D marginal, red line at the median.
Lower triangle: 2D density. Upper triangle: Pearson correlation. Four chains pooled,
16,000 post-warmup draws per group.

Medians for the same fit are in ../data/posteriors_invivo_withDIP.csv.

## Supporting robustness figures (no-DIP comparison)

Referenced in the manuscript Supporting Information subsection "Independent test: a
no-DIP model with both positive and negative IFN feedback".

- `invivo_noDIP_group_specific_cumF50_fits.pdf` — no-DIP model fit to in vivo V and F,
  with V0 and cumF50 allowed to vary by group and all other kinetic parameters shared.
  Posterior median and 50/90% credible intervals; points are individual mice, open
  triangles are below-detection observations.
- `invivo_withDIP_vs_noDIP_crossfit.png` — posterior-predictive overlay of the with-DIP
  and no-DIP fits (both with cumulative-IFN feedback) on the in vivo BALB/c mouse data.
  Rows: low/med/high DVG inoculum; columns: V and F. Coloured = with-DIP, black = no-DIP.
