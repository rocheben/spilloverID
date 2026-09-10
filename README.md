# spilloverID: Reservoir Identification for Zoonotic Spillover Without Genetic Data

An R implementation of a composite epidemiological plausibility score (Π) that
ranks candidate wildlife hosts as the source of human infection using only
serological and exposure data — no genetic data required. The score integrates
14 independent epidemiological signatures across spatial, individual,
multi-level, demographic, temporal, mechanistic and specificity dimensions,
together with three maintenance signatures (C12–C14) that separate a maintenance
reservoir from an amplifying or bridge host.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rocheben/spilloverID")
```

or, from the deposited archive:

```sh
R CMD INSTALL package/spilloverID
```

## Quick start

```r
library(spilloverID)

# 1. Build a spillover_data object from your field data
data <- spillover_data(
  human        = my_human_survey,    # individual-level data with exposure columns
  animal       = my_animal_summary,  # site x host prevalence
  animal_indiv = my_animal_indiv,    # optional, for the catalytic models
  site         = my_site_covariates,
  n_hosts      = 4,
  host_names   = c("Oryzomys", "Peromyscus", "Heteromys", "Mus")
)

# 2. Compute the composite plausibility score
result <- compute_pi(data, bootstrap = 100)
print(result)
```

`compute_pi()` returns the score per host, the per-criterion matrix, the weights
actually applied after redistribution, the ranks, and the bootstrap matrix when
one was requested.

## What to check before acting on a ranking

- **Availability.** `result$weights_used` shows which criteria were evaluable. A
  criterion returns `NA` either because its input data are absent or because its
  fit did not converge, and its weight is redistributed over the rest. The scale
  of Π is not comparable between two runs with different availability patterns.
- **Margin, not rank.** The gap between the first- and second-ranked host, with
  the percentile interval of that difference across the bootstrap replicates, is
  what tells you whether the ranking is worth acting on.
- **C12 against the overall ranking.** A divergence between the two is the
  signature of a bridge-host or vector-borne configuration.
- **The candidate set.** Π ranks the hosts you sampled against one another. It
  provides no test that the set contains a reservoir at all.

## Simulation engine

```r
# Score a dataset whose reservoir is known by construction
sim <- simulate_spillover(scenario_params("S1_canonical"), seed = 42)
attr(sim, "truth")$k_star
compute_pi(sim)
```

## Vignettes

- `vignette("getting-started", package = "spilloverID")` — the field-data workflow
- `vignette("benchmark", package = "spilloverID")` — the simulation benchmark
- `vignette("case-studies", package = "spilloverID")` — the three reconstructions

## Citation

Roche, B., Bento, A., Drexler, F. & Arnal, A. (2026). Prediction of wildlife
reservoirs of zoonotic pathogens from observational data: a validated composite
scoring framework.

## License

MIT
