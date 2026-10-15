# ============================================================================
# Comparability of the Pi ~ 0.55 interpretation threshold.
#
# The threshold is calibrated on the simulated benchmark, where the negative
# control S11 has 13 ACTIVE criteria (only C11 returns NA). It is then applied
# to the three real cases, where C8 AND C11 return NA -- 12 active criteria,
# with weight redistributed over a smaller set. In addition C9 returns exactly
# 0 for every host in all three real cases, so its effective weight (0.065) is
# dead: it depresses every Pi by a constant factor without discriminating.
#
# This script recomputes the S11 max-Pi null distribution under four
# criterion-availability regimes, matched replicate-by-replicate (same seeds,
# same simulated datasets), and reports the threshold in each regime that
# preserves the false-positive rate of Pi >= 0.55 in the benchmark regime.
#   bench13  : the simulated benchmark (C11 NA)            -> 13 active
#   real12   : C8 forced NA as well                        -> 12 active
#   real12_c9dead : C8 NA and C9 forced to 0 for all hosts -> 12 active, 1 dead
#   real11   : C8 and C9 both forced NA                    -> 11 active
# ============================================================================
source("../_load_cases.R")
dir.create("../../../data/suppl_S5_uncertainty", showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(parallel))

N_REP  <- 300
NCORES <- { n <- suppressWarnings(detectCores()); if (!is.finite(n)) n <- as.integer(Sys.getenv("NPROC", "2")); max(1L, as.integer(Sys.getenv("S5_N_CORE", as.character(n)))) }
PUB_THRESHOLD <- 0.55

# compute_pi with selected criteria forced to NA (drop = weight redistributed)
# or forced to 0 (dead = weight retained but non-discriminating).
compute_pi_masked <- function(data, drop = character(0), zero = character(0),
                              weights = NULL, crit_mat = NULL) {
  if (is.null(weights)) weights <- default_weights()
  weights <- weights[names(weights) %in% paste0("C", 1:14)]
  weights <- weights / sum(weights)
  crit_funs <- list(
    C1 = criterion_C1_spatial, C2 = criterion_C2_reservoir_signatures,
    C3 = criterion_C3_individual, C4 = criterion_C4_dose_response,
    C5 = criterion_C5_cryptic_fraction, C6 = criterion_C6_interaction,
    C7 = criterion_C7_foi_coherence, C8 = criterion_C8_temporal_precedence,
    C9 = criterion_C9_mediation, C10 = criterion_C10_mechanistic,
    C11 = criterion_C11_specificity, C12 = criterion_C12_inter_epidemic,
    C13 = criterion_C13_endemic_age, C14 = criterion_C14_spatial_independence)
  K <- data$n_hosts
  # The criterion matrix does not depend on the regime: the regimes only mask it
  # afterwards. It is therefore computed once per replicate and passed in, which
  # makes the four regimes cost one criterion pass instead of four.
  if (is.null(crit_mat)) {
    crit_mat <- matrix(NA_real_, nrow = K, ncol = length(crit_funs),
                       dimnames = list(data$host_names, names(crit_funs)))
    for (cn in names(crit_funs))
      crit_mat[, cn] <- tryCatch(crit_funs[[cn]](data),
                                 error = function(e) rep(NA_real_, K))
  }
  crit_mat[, colnames(crit_mat) %in% drop] <- NA_real_
  crit_mat[, colnames(crit_mat) %in% zero] <- 0
  na_crit <- apply(crit_mat, 2, function(x) all(is.na(x)))
  w_active <- weights
  w_active[names(w_active) %in% names(na_crit[na_crit])] <- 0
  if (sum(w_active) > 0) w_active <- w_active / sum(w_active)
  crit_safe <- crit_mat; crit_safe[is.na(crit_safe)] <- 0
  pi_vec <- as.numeric(crit_safe %*% w_active)
  names(pi_vec) <- data$host_names
  list(pi = pi_vec, criteria = crit_mat, weights_used = w_active,
       n_active = sum(w_active > 0))
}

# The regimes differ in which criteria are available. `reconstruction11`
# reproduces the availability pattern of the three outbreak reconstructions,
# where C6 does not converge, C8 has no time series and C11 has no negative
# control, leaving eleven active criteria with C9 present but returning zero.
REGIMES <- list(
  bench13          = list(drop = character(0),   zero = character(0)),
  real12           = list(drop = "C8",           zero = character(0)),
  real12_c9dead    = list(drop = "C8",           zero = "C9"),
  real11           = list(drop = c("C8", "C9"),  zero = character(0)),
  reconstruction11 = list(drop = c("C6", "C8"),  zero = "C9")
)

# Seeds identical to run_benchmark.R for S11 so the null datasets match the
# published negative control replicate for replicate.
scenarios <- c("S1_canonical","S2_bridge","S3_cryptic","S4_confounded",
               "S5_codistributed","S6_strong_seasonal","S7_low_power",
               "S8_measurement_error","S9_weak_signal","S10_short_immunity",
               "S11_no_spillover")
seed_for <- function(rep, scen) 20260524 + rep * 17 + match(scen, scenarios)

one_rep <- function(rep, scen) {
  d <- tryCatch(simulate_spillover(scenario_params(scen), seed = seed_for(rep, scen)),
                error = function(e) NULL)
  if (is.null(d)) return(NULL)
  k_star <- attr(d, "truth")$k_star
  base <- tryCatch(suppressWarnings(compute_pi_masked(d)), error = function(e) NULL)
  if (is.null(base)) return(NULL)
  out <- list()
  for (rg in names(REGIMES)) {
    r <- tryCatch(suppressWarnings(compute_pi_masked(d, drop = REGIMES[[rg]]$drop,
                                                     zero = REGIMES[[rg]]$zero,
                                                     crit_mat = base$criteria)),
                  error = function(e) NULL)
    if (is.null(r)) next
    srt <- sort(r$pi, decreasing = TRUE)
    out[[rg]] <- data.frame(
      scenario = scen, rep = rep, regime = rg, n_active = r$n_active,
      max_pi = as.numeric(srt[1]), second_pi = as.numeric(srt[2]),
      margin_top12 = as.numeric(srt[1] - srt[2]),
      pi_true = if (is.na(k_star)) NA_real_ else as.numeric(r$pi[k_star]),
      rank_true = if (is.na(k_star)) NA_integer_ else as.integer(rank(-r$pi, ties.method = "min")[k_star]),
      margin_true = if (is.na(k_star)) NA_real_ else
        as.numeric(r$pi[k_star] - max(r$pi[-k_star])),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# ---- negative control S11 under all four regimes ----------------------------
cat("S11 negative control,", N_REP, "replicates x", length(REGIMES), "regimes\n"); flush.console()
t0 <- Sys.time()
s11 <- mclapply(seq_len(N_REP), one_rep, scen = "S11_no_spillover", mc.cores = NCORES)
s11 <- do.call(rbind, s11)
cat("  elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---- positive scenarios under bench13 vs real12 (scale shift on true hosts) --
# A subset is enough to characterise the scale shift; use the two extremes of
# performance (S1 best, S8 worst) plus the canonical mid scenarios.
POS <- c("S1_canonical", "S7_low_power", "S8_measurement_error")
N_POS  <- as.integer(Sys.getenv("S5_N_POS", "50"))
pos <- list()
for (sc in POS) {
  cat("positive scenario", sc, "\n"); flush.console()
  p <- mclapply(seq_len(N_POS), one_rep, scen = sc, mc.cores = NCORES)
  pos[[sc]] <- do.call(rbind, p)
}
pos <- do.call(rbind, pos)

allrep <- rbind(s11, pos)
write.csv(allrep, "../../../data/suppl_S5_uncertainty/regime_per_rep.csv", row.names = FALSE)

# ---- threshold recalibration ------------------------------------------------
# Two equivalent-threshold definitions, both anchored on the benchmark regime:
#  (a) FPR-matched: the quantile of the regime's S11 max-Pi null distribution at
#      the same tail probability that 0.55 occupies in the bench13 null.
#  (b) quantile-mapped: q_regime(F_bench13(0.55)) -- identical construction,
#      reported for both 0.55 and the 0.50 / 0.60 sensitivity thresholds.
nulls <- split(s11$max_pi, s11$regime)
bench_null <- nulls[["bench13"]]

equiv_threshold <- function(thr, bench, target) {
  p <- mean(bench >= thr)                 # FPR in the calibration regime
  if (p <= 0) return(max(target))         # threshold beyond the observed null
  as.numeric(quantile(target, 1 - p, type = 7))
}

trows <- list()
for (rg in names(nulls)) {
  tg <- nulls[[rg]]
  for (thr in c(0.50, 0.55, 0.60)) {
    trows[[length(trows) + 1]] <- data.frame(
      regime = rg,
      n_active = median(s11$n_active[s11$regime == rg]),
      published_threshold = thr,
      fpr_bench13 = mean(bench_null >= thr),
      fpr_this_regime_at_published = mean(tg >= thr),
      equivalent_threshold = equiv_threshold(thr, bench_null, tg),
      null_median = median(tg), null_q95 = quantile(tg, 0.95),
      null_q975 = quantile(tg, 0.975), null_max = max(tg),
      stringsAsFactors = FALSE)
  }
}
thresholds <- do.call(rbind, trows)
write.csv(thresholds, "../../../data/suppl_S5_uncertainty/threshold_recalibration.csv", row.names = FALSE)

cat("\n--- S11 null distribution by regime ---\n")
print(aggregate(max_pi ~ regime + n_active, data = s11,
                FUN = function(x) c(med = median(x), q95 = quantile(x, .95),
                                    fp55 = mean(x >= .55))), digits = 3)
cat("\n--- threshold recalibration ---\n")
print(thresholds[thresholds$published_threshold == 0.55, ], row.names = FALSE, digits = 3)
cat("\n--- positive scenarios: pi_true shift bench13 -> real12 ---\n")
print(aggregate(pi_true ~ scenario + regime, data = pos, FUN = median), digits = 3)
cat("\n--- top1 rate by scenario and regime ---\n")
pos$top1 <- as.integer(pos$rank_true == 1)
print(aggregate(top1 ~ scenario + regime, data = pos, FUN = mean), digits = 3)
saveRDS(list(s11 = s11, pos = pos, thresholds = thresholds), "../../../data/suppl_S5_uncertainty/regime_raw.rds")
cat("\nDONE 03\n")
