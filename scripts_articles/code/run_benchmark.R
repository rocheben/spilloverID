#' ============================================================================
#' Mechanistic benchmark for spilloverID: 11 scenarios x N_REP replicates.
#'
#' A single pass over the simulated datasets records everything the manuscript
#' and the supplementary sections need, so that no scenario is ever simulated
#' twice:
#'   - the composite score, its rank and the oracle discrimination margin
#'     (Pi of the true reservoir minus the highest Pi among the others);
#'   - the OBSERVABLE margin (highest Pi minus second-highest), which is what a
#'     field team can actually compute, used by the interpretation rules and by
#'     Supplementary Section S5;
#'   - the full per-criterion matrix, used by the composite-versus-single-
#'     criterion comparison (Supplementary Figure 4);
#'   - the reference rankings of Supplementary Section S6.
#'
#' Writes ../data/benchmark_results.json      (summary + per-replicate)
#'        ../data/benchmark_criteria.csv      (per-replicate criterion scores)
#'        ../data/benchmark_baselines.csv     (per-replicate reference rankings)
#' Checkpointed per scenario: restarting resumes where it stopped.
#' Runtime: several hours on two cores at N_REP = 300.
#' ============================================================================

suppressPackageStartupMessages({
  library(spilloverID)
  library(jsonlite)
  library(parallel)
})

source("baselines.R")

N_REP  <- as.integer(Sys.getenv("BENCH_N_REP", "300"))
.dc    <- suppressWarnings(parallel::detectCores())
if (!is.finite(.dc)) .dc <- 1L
N_CORE <- max(1L, as.integer(Sys.getenv("BENCH_N_CORE", as.character(.dc))))

scenarios <- c("S1_canonical", "S2_bridge", "S3_cryptic", "S4_confounded",
               "S5_codistributed", "S6_strong_seasonal", "S7_low_power",
               "S8_measurement_error", "S9_weak_signal", "S10_short_immunity",
               "S11_no_spillover")

scenario_labels <- c(
  "S1_canonical"         = "S1 canonical",
  "S2_bridge"            = "S2 bridge host",
  "S3_cryptic"           = "S3 cryptic",
  "S4_confounded"        = "S4 confounded",
  "S5_codistributed"     = "S5 co-distributed",
  "S6_strong_seasonal"   = "S6 seasonal forcing",
  "S7_low_power"         = "S7 low power",
  "S8_measurement_error" = "S8 measurement error",
  "S9_weak_signal"       = "S9 weak signal",
  "S10_short_immunity"   = "S10 short immunity",
  "S11_no_spillover"     = "S11 no spillover (negative control)")

ckpt_dir <- "../data/benchmark_checkpoints"
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

#' One replicate: simulate, score, and record every downstream quantity.
one_rep <- function(scen, rep) {
  seed <- 20260524 + rep * 17 + match(scen, scenarios)
  d <- tryCatch(simulate_spillover(scenario_params(scen), seed = seed),
                error = function(e) NULL)
  if (is.null(d)) return(NULL)
  res <- tryCatch(suppressWarnings(compute_pi(d, bootstrap = 0)),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)

  k_star  <- attr(d, "truth")$k_star
  pi_vec  <- res$pi
  srt     <- sort(pi_vec, decreasing = TRUE)
  obs_margin <- as.numeric(srt[1] - srt[2])   # observable: top1 - top2

  base <- data.frame(
    scenario = scen, rep = rep,
    k_star = if (is.na(k_star)) NA_integer_ else as.integer(k_star),
    pi_true = if (is.na(k_star)) NA_real_ else as.numeric(pi_vec[k_star]),
    pi_max_other = if (is.na(k_star)) as.numeric(srt[1])
                   else as.numeric(max(pi_vec[-k_star])),
    max_pi = as.numeric(srt[1]),
    pi_range = as.numeric(srt[1] - srt[length(srt)]),
    obs_margin = obs_margin,
    pred_top1 = as.integer(which.max(pi_vec)),
    stringsAsFactors = FALSE)
  base$margin     <- if (is.na(k_star)) NA_real_ else base$pi_true - base$pi_max_other
  base$rank_true  <- if (is.na(k_star)) NA_integer_ else as.integer(res$rank[k_star])
  base$top1       <- if (is.na(k_star)) NA_integer_ else as.integer(base$rank_true == 1)
  base$top2       <- if (is.na(k_star)) NA_integer_ else as.integer(base$rank_true <= 2)

  crit <- as.data.frame(res$criteria)
  crit$scenario <- scen; crit$rep <- rep; crit$host <- rownames(res$criteria)
  crit$is_true <- if (is.na(k_star)) NA else seq_len(nrow(crit)) == k_star

  scores <- all_baseline_scores(d, pi_vec, res$criteria)
  bl <- do.call(rbind, lapply(names(scores), function(nm) {
    r <- if (is.na(k_star)) list(rank = NA_integer_, top1_strict = NA_integer_,
                                 top1_exp = NA_real_, top2 = NA_integer_,
                                 n_tied_top = NA_integer_)
         else rank_true_host(scores[[nm]], k_star)
    data.frame(scenario = scen, rep = rep, baseline = nm,
               rank = r$rank, top1_strict = r$top1_strict,
               top1_exp = r$top1_exp, top2 = r$top2,
               n_tied_top = r$n_tied_top,
               top_gap = top_gap(scores[[nm]]),
               stringsAsFactors = FALSE)
  }))

  list(base = base, crit = crit, bl = bl)
}

for (scen in scenarios) {
  ck <- file.path(ckpt_dir, paste0(scen, ".rds"))
  if (file.exists(ck)) { cat("[", scen, "] cached\n", sep = ""); next }
  cat("[", scen, "] ", N_REP, " replicates on ", N_CORE, " core(s) ... ",
      sep = ""); flush.console()
  t0 <- Sys.time()
  got <- mclapply(seq_len(N_REP), function(rep) one_rep(scen, rep),
                  mc.cores = N_CORE)
  got <- Filter(function(x) !is.null(x) && !inherits(x, "try-error"), got)
  saveRDS(got, ck)
  cat(length(got), "ok in",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}

all_got <- unlist(lapply(scenarios,
  function(s) readRDS(file.path(ckpt_dir, paste0(s, ".rds")))), recursive = FALSE)

per_rep <- do.call(rbind, lapply(all_got, `[[`, "base"))
crit_df <- do.call(rbind, lapply(all_got, `[[`, "crit"))
bl_df   <- do.call(rbind, lapply(all_got, `[[`, "bl"))

per_rep$`_row` <- paste0(per_rep$scenario, ".h", ifelse(is.na(per_rep$k_star), 1, per_rep$k_star))

summary_rows <- lapply(scenarios, function(scen) {
  s <- per_rep[per_rep$scenario == scen, ]
  neg <- scen == "S11_no_spillover"
  out <- list(scenario = scen, label = unname(scenario_labels[scen]),
              n_rep = nrow(s))
  if (!neg) {
    out$top1_rate <- round(100 * mean(s$top1, na.rm = TRUE), 1)
    out$top2_rate <- round(100 * mean(s$top2, na.rm = TRUE), 1)
    out$mean_margin <- round(mean(s$margin, na.rm = TRUE), 3)
    out$sd_margin <- round(stats::sd(s$margin, na.rm = TRUE), 3)
    out$pct_margin_negative <- round(100 * mean(s$margin < 0, na.rm = TRUE), 1)
    out$median_pi_true <- round(stats::median(s$pi_true, na.rm = TRUE), 3)
  }
  out$median_max_pi <- round(stats::median(s$max_pi, na.rm = TRUE), 3)
  out$median_pi_range <- round(stats::median(s$pi_range, na.rm = TRUE), 3)
  out$median_obs_margin <- round(stats::median(s$obs_margin, na.rm = TRUE), 3)
  out$`_row` <- scen
  out
})

write_json(list(n_rep_per_scenario = N_REP,
                summary = summary_rows, per_rep = per_rep),
           "../data/benchmark_results.json",
           auto_unbox = TRUE, pretty = TRUE, digits = 4, na = "null")
write.csv(crit_df, "../data/benchmark_criteria.csv", row.names = FALSE)
write.csv(bl_df,   "../data/benchmark_baselines.csv", row.names = FALSE)

cat("\n== DONE ==\n")
cat("Wrote benchmark_results.json, benchmark_criteria.csv, benchmark_baselines.csv\n")
