#' ============================================================================
#' Observable versus oracle discrimination margin.
#'
#' Derived from the per-replicate benchmark output rather than from a second
#' simulation pass: run_benchmark.R already records, for every replicate, the
#' oracle margin (true reservoir minus best competitor, computable only when the
#' truth is known) and the observable margin (highest minus second-highest
#' score, which is what a field team can compute). This script reshapes those
#' columns into the file consumed by 05_calibration.py.
#'
#' Requires ../../../data/benchmark_results.json
#' Writes   ../../../data/suppl_S5_uncertainty/benchmark_recomputed_margins.csv
#' Runtime: seconds.
#' ============================================================================
suppressPackageStartupMessages(library(jsonlite))

OUT <- "../../../data/suppl_S5_uncertainty"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

b <- fromJSON("../../../data/benchmark_results.json")$per_rep
bench <- data.frame(
  scenario  = b$scenario,
  rep       = b$rep,
  n_active  = NA_integer_,
  k_star    = b$k_star,
  max_pi    = b$max_pi,
  second_pi = b$max_pi - b$obs_margin,
  margin_obs = b$obs_margin,
  pi_true   = b$pi_true,
  rank_true = b$rank_true,
  margin_oracle = b$margin,
  top1      = b$top1,
  stringsAsFactors = FALSE)

write.csv(bench, file.path(OUT, "benchmark_recomputed_margins.csv"), row.names = FALSE)
cat("Wrote benchmark_recomputed_margins.csv (", nrow(bench), " replicates )\n", sep = "")
