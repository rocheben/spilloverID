#' ============================================================================
#' Supplementary Section S6 -- the composite score against reference rankings.
#'
#' Eight rankings computable from exactly the same input data as Pi are scored
#' on the simulated benchmark and on the three outbreak reconstructions:
#'   sero  pooled animal seroprevalence          pcr   pooled animal PCR prevalence
#'   or    crude declared-exposure odds ratio    C12   inter-epidemic persistence alone
#'   C7    catalytic force-of-infection alone    C14   spatial independence alone
#'   C5    declared-exposure fraction alone      taxo  a-priori taxonomic-order prior
#'
#' The taxonomic prior is an explicitly labelled caricature of the only element
#' of trait- and phylogeny-based predictors that transfers to a within-outbreak
#' setting; it assigns the same score to every species of an order, so it is
#' structurally tied on the simulated candidate sets and is reported on the three
#' reconstructions only.
#'
#' The simulated columns are read from benchmark_baselines.csv, written by
#' run_benchmark.R over the same replicates as the benchmark itself.
#'
#' Requires ../../data/benchmark_baselines.csv
#' Writes   ../../data/suppl_S6_baselines/
#' Runtime: ~2 min (the case pass).
#' ============================================================================

source("../_load_cases.R")
source("../../baselines.R")
suppressPackageStartupMessages(library(jsonlite))

OUT <- "../../../data/suppl_S6_baselines"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- simulated side, pooled over the ten positive scenarios -----------------
bl <- utils::read.csv("../../../data/benchmark_baselines.csv", stringsAsFactors = FALSE)
pos <- bl[bl$scenario != "S11_no_spillover", ]
neg <- bl[bl$scenario == "S11_no_spillover", ]

pooled <- do.call(rbind, lapply(BASELINE_IDS, function(b) {
  s <- pos[pos$baseline == b, ]
  n <- neg[neg$baseline == b, ]
  data.frame(baseline = b,
             n_rep_per_scenario = length(unique(s$rep)),
             n_eff = sum(!is.na(s$rank)),
             sim_top1_pct = round(100 * mean(s$top1_strict, na.rm = TRUE), 1),
             sim_top1_pct_tiebroken = round(100 * mean(s$top1_exp, na.rm = TRUE), 1),
             sim_top2_pct = round(100 * mean(s$top2, na.rm = TRUE), 1),
             S11_median_top_gap = round(stats::median(n$top_gap, na.rm = TRUE), 3),
             `S11_pct_gap_gt_0.5` = round(100 * mean(n$top_gap > 0.5, na.rm = TRUE), 1),
             check.names = FALSE)
}))

per_scen <- do.call(rbind, lapply(split(pos, list(pos$scenario, pos$baseline)),
  function(s) if (nrow(s) == 0) NULL else data.frame(
    scenario = s$scenario[1], baseline = s$baseline[1], n_eff = nrow(s),
    top1_strict = round(100 * mean(s$top1_strict, na.rm = TRUE), 1),
    top1_exp = round(100 * mean(s$top1_exp, na.rm = TRUE), 1),
    top2 = round(100 * mean(s$top2, na.rm = TRUE), 1))))
rownames(per_scen) <- NULL

# ---- the three reconstructions ---------------------------------------------
cases <- build_cases()
case_rows <- list()
for (k in names(cases)) {
  d   <- cases[[k]]
  res <- suppressWarnings(compute_pi(d, bootstrap = 0))
  kt  <- which(names(res$pi) == CASE_TRUE[[k]])
  sc  <- all_baseline_scores(d, res$pi, res$criteria)
  for (b in BASELINE_IDS) {
    r <- rank_true_host(sc[[b]], kt)
    case_rows[[length(case_rows) + 1]] <- data.frame(
      case = k, case_label = unname(CASE_LAB[[k]]), baseline = b,
      true_host = CASE_TRUE[[k]], rank_true = r$rank,
      top1_strict = r$top1_strict, top2 = r$top2, n_tied_top = r$n_tied_top)
  }
}
case_ranks <- do.call(rbind, case_rows)

agg <- aggregate(cbind(top1_strict, top2) ~ baseline, data = case_ranks,
                 FUN = function(x) sum(x, na.rm = TRUE), na.action = stats::na.pass)
names(agg) <- c("baseline", "real_cases_top1_of3", "real_cases_top2_of3")
rk <- reshape(case_ranks[, c("case", "baseline", "rank_true")],
              idvar = "baseline", timevar = "case", direction = "wide")
names(rk) <- sub("^rank_true\\.", "rank_", names(rk))

tab <- Reduce(function(a, b) merge(a, b, by = "baseline"), list(pooled, agg, rk))
tab$baseline_label <- unname(BASELINE_LABELS[tab$baseline])
tab <- tab[match(BASELINE_IDS, tab$baseline), ]
tab <- tab[, c("baseline", "baseline_label", "n_rep_per_scenario",
               "sim_top1_pct", "sim_top1_pct_tiebroken", "sim_top2_pct", "n_eff",
               "S11_median_top_gap", "S11_pct_gap_gt_0.5",
               "rank_sn", "rank_nv", "rank_rv",
               "real_cases_top1_of3", "real_cases_top2_of3")]

write.csv(tab,        file.path(OUT, "pi_vs_baselines_table.csv"), row.names = FALSE)
write.csv(per_scen,   file.path(OUT, "bench_per_scenario.csv"),    row.names = FALSE)
write.csv(case_ranks, file.path(OUT, "case_baseline_ranks.csv"),   row.names = FALSE)
print(tab[, c("baseline", "sim_top1_pct", "sim_top2_pct",
              "rank_sn", "rank_nv", "rank_rv", "real_cases_top1_of3")],
      row.names = FALSE)
cat("S6 DONE\n")
