#' ============================================================================
#' Supplementary Section S7 -- taxonomic resolution, multiple reservoirs and
#' inventory completeness.
#'
#' Extends the mechanistic benchmark with scenarios that relax the assumption of
#' a single, unambiguously delimited reservoir species: two co-located congeners
#' of which only one seeds human infection (with a cross-reactive and with a
#' species-resolving assay), two hosts that genuinely co-maintain the pathogen
#' (symmetric and asymmetric spillover), and a candidate set from which the true
#' reservoir is absent. Two reference arms reproduce the single-reservoir and
#' no-spillover configurations under the same five-host geometry.
#'
#' Writes ../../data/suppl_S7_taxonomic/. Checkpointed per scenario.
#
# Usage:  Rscript run_benchmark_ext.R [N_REP] [N_CORES]
# Seeds:  seed(scenario, rep) = 20260524 + 17 * rep + 1000 * scenario_index
#         (offset by 1000 per scenario so no seed collides with the
#          published benchmark, which uses 20260524 + 17*rep + index_1_11)
# =====================================================================
suppressPackageStartupMessages({ library(spilloverID); library(jsonlite); library(parallel) })

args    <- commandArgs(trailingOnly = TRUE)
.dc     <- suppressWarnings(detectCores()); if (!is.finite(.dc)) .dc <- 1L
N_REP   <- if (length(args) >= 1) as.integer(args[1]) else
             as.integer(Sys.getenv("EXT_N_REP", "150"))
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else
             as.integer(Sys.getenv("EXT_N_CORE", as.character(.dc)))

source("simulate_ext.R")
source("eval_ext.R")

OUT <- "../../../data/suppl_S7_taxonomic"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SCENARIOS <- c("S12_congeners", "S12b_congeners_resolved",
               "S13_dual_symmetric", "S13b_dual_asymmetric",
               "S14_hidden_reservoir",
               "S1c_ref_single", "S11c_ref_null")

SCEN_LABELS <- c(
  S12_congeners            = "S12 congeners (cross-reactive assay)",
  S12b_congeners_resolved  = "S12b congeners (species-resolving assay)",
  S13_dual_symmetric       = "S13 dual reservoir (symmetric)",
  S13b_dual_asymmetric     = "S13b dual reservoir (asymmetric)",
  S14_hidden_reservoir     = "S14 incomplete inventory (reservoir unsampled)",
  S1c_ref_single           = "S1c reference: single reservoir",
  S11c_ref_null            = "S11c reference: no spillover"
)

seed_for <- function(scen, rep) 20260524L + 17L * rep + 1000L * match(scen, SCENARIOS)

# ---- one replicate -------------------------------------------------
run_one <- function(scen, rep) {
  sd_use <- seed_for(scen, rep)
  d <- tryCatch(simulate_multihost(scenario_params_ext(scen), seed = sd_use),
                error = function(e) NULL)
  if (is.null(d)) return(NULL)
  tr <- attr(d, "truth")
  res <- tryCatch(suppressWarnings(compute_pi(d, bootstrap = 0)),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  ev <- eval_replicate(res$pi, tr, scen, rep, sd_use)
  # per-criterion contrast between the two focal hosts, where meaningful
  ctr <- NULL
  if (scen %in% c("S12_congeners","S12b_congeners_resolved")) {
    ctr <- criteria_contrast(res$criteria, 1, 2, scen, rep)   # true vs congener
  } else if (scen %in% c("S13_dual_symmetric","S13b_dual_asymmetric")) {
    ctr <- criteria_contrast(res$criteria, 1, 2, scen, rep)   # co-reservoir 1 vs 2
  }
  list(ev = ev, ctr = ctr,
       wts = res$weights_used,
       na_crit = paste(colnames(res$criteria)[
         apply(res$criteria, 2, function(x) all(is.na(x)))], collapse = ","))
}

# ---- loop over scenarios -------------------------------------------
all_ev <- list(); all_ctr <- list(); meta <- list()
t_start <- Sys.time()

for (scen in SCENARIOS) {
  ck <- file.path(OUT, paste0("ckpt_", scen, ".rds"))
  if (file.exists(ck)) {
    got <- readRDS(ck)
    if (nrow(got$ev) >= N_REP) {
      cat("[", scen, "] checkpoint complete (", nrow(got$ev), " reps)\n", sep = "")
      all_ev[[scen]] <- got$ev; all_ctr[[scen]] <- got$ctr; meta[[scen]] <- got$meta
      next
    }
  }
  cat("[", scen, "] running ", N_REP, " reps on ", N_CORES, " cores ... ", sep = "")
  t0 <- Sys.time()
  out <- mclapply(seq_len(N_REP), function(r) run_one(scen, r),
                  mc.cores = N_CORES, mc.preschedule = TRUE)
  out <- Filter(function(z) !is.null(z) && !inherits(z, "try-error"), out)
  ev  <- do.call(rbind, lapply(out, `[[`, "ev"))
  ctr_list <- Filter(Negate(is.null), lapply(out, `[[`, "ctr"))
  ctr <- if (length(ctr_list)) do.call(rbind, ctr_list) else NULL
  m <- list(n_requested = N_REP, n_completed = nrow(ev),
            na_criteria = names(sort(table(sapply(out, `[[`, "na_crit")),
                                     decreasing = TRUE))[1],
            weights_used = out[[1]]$wts,
            seed_formula = "20260524 + 17*rep + 1000*scenario_index",
            seed_first = seed_for(scen, 1), seed_last = seed_for(scen, N_REP),
            elapsed_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(list(ev = ev, ctr = ctr, meta = m), ck)
  all_ev[[scen]] <- ev; all_ctr[[scen]] <- ctr; meta[[scen]] <- m
  cat(nrow(ev), " ok in ", round(m$elapsed_s), "s\n", sep = "")
}

# Scenarios differ in K (4 or 5 observed hosts), so the per-replicate
# frames carry different pi_hK columns. Align on the union, filling NA.
rbind_fill <- function(lst) {
  lst <- Filter(function(x) !is.null(x) && nrow(x) > 0, lst)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[, cols, drop = FALSE]
  }))
}
ev_all  <- rbind_fill(all_ev)
ctr_all <- rbind_fill(all_ctr)
if (!is.null(ctr_all) && nrow(ctr_all) == 0) ctr_all <- NULL

# ---- summarise ------------------------------------------------------
pct <- function(x) if (all(is.na(x))) NA_real_ else round(mean(x, na.rm = TRUE) * 100, 1)
r3  <- function(x) if (all(is.na(x))) NA_real_ else round(mean(x, na.rm = TRUE), 3)
m3  <- function(x) if (all(is.na(x))) NA_real_ else round(stats::median(x, na.rm = TRUE), 3)

summ <- do.call(rbind, lapply(SCENARIOS, function(scen) {
  d <- ev_all[ev_all$scenario == scen, ]
  data.frame(
    scenario         = scen,
    label            = unname(SCEN_LABELS[scen]),
    n_rep            = nrow(d),
    n_hosts          = d$n_hosts[1],
    n_true_reservoirs= d$n_true_reservoirs[1],
    top1_species     = pct(d$top1_species),
    top1_genus       = pct(d$top1_genus),
    both_top2        = pct(d$both_top2),
    any_top1         = pct(d$any_top1),
    mean_margin      = r3(d$margin),
    mean_margin_set  = r3(d$margin_set),
    mean_gap_within  = r3(d$pi_gap_within_true),
    median_pi_true   = m3(d$pi_true),
    median_max_pi    = m3(d$max_pi),
    median_pi_range  = m3(d$pi_range),
    pct_max_pi_gt_055= pct(d$exceeds_055),
    stringsAsFactors = FALSE
  )
}))

cat("\n================ EXTENDED BENCHMARK ================\n")
print(summ[, c("scenario","n_rep","top1_species","top1_genus","both_top2",
               "mean_margin","median_max_pi","pct_max_pi_gt_055")], row.names = FALSE)

# ---- S12: does species resolution depend on ecological similarity? --
sim_fit <- list()
for (scen in c("S12_congeners","S12b_congeners_resolved")) {
  d <- ev_all[ev_all$scenario == scen, ]
  if (nrow(d) > 10 && length(unique(d$top1_species)) > 1) {
    f <- glm(top1_species ~ rho_congener, data = d, family = binomial())
    cf <- summary(f)$coefficients
    # predicted P(exact species) at rho = 0.70 and 0.90
    pr <- predict(f, newdata = data.frame(rho_congener = c(0.70, 0.80, 0.90)),
                  type = "response")
    sim_fit[[scen]] <- list(
      slope = unname(cf["rho_congener","Estimate"]),
      se    = unname(cf["rho_congener","Std. Error"]),
      p     = unname(cf["rho_congener","Pr(>|z|)"]),
      p_exact_rho070 = unname(pr[1]), p_exact_rho080 = unname(pr[2]),
      p_exact_rho090 = unname(pr[3]))
  } else {
    sim_fit[[scen]] <- list(slope = NA, se = NA, p = NA,
      p_exact_rho070 = NA, p_exact_rho080 = NA, p_exact_rho090 = NA)
  }
}

# ---- S13b: does the score track the true share, or arbitrate? -------
d13b <- ev_all[ev_all$scenario == "S13b_dual_asymmetric", ]
share_fit <- if (nrow(d13b) > 10) {
  f <- lm(pi_gap_within_true ~ share_1, data = d13b)
  list(slope = unname(coef(f)["share_1"]),
       p = summary(f)$coefficients["share_1","Pr(>|t|)"],
       cor_gap_share = unname(cor(d13b$pi_gap_within_true, d13b$share_1)),
       correct_dominant = mean(d13b$pi_h1 > d13b$pi_h2))
} else list(slope = NA, p = NA, cor_gap_share = NA, correct_dominant = NA)

# ---- S14 vs S11c: specificity of the score under incomplete inventory
d14  <- ev_all[ev_all$scenario == "S14_hidden_reservoir", ]
d11c <- ev_all[ev_all$scenario == "S11c_ref_null", ]
d1c  <- ev_all[ev_all$scenario == "S1c_ref_single", ]
spec <- list(
  S14_median_max_pi  = m3(d14$max_pi),
  S14_q025 = round(quantile(d14$max_pi, 0.025), 3),
  S14_q975 = round(quantile(d14$max_pi, 0.975), 3),
  S11c_median_max_pi = m3(d11c$max_pi),
  S11c_q975 = round(quantile(d11c$max_pi, 0.975), 3),
  S1c_median_max_pi  = m3(d1c$max_pi),
  S14_pct_gt_055     = pct(d14$exceeds_055),
  S11c_pct_gt_055    = pct(d11c$exceeds_055),
  S1c_pct_gt_055     = pct(d1c$exceeds_055),
  wilcox_S14_vs_S11c_p = tryCatch(
    wilcox.test(d14$max_pi, d11c$max_pi)$p.value, error = function(e) NA),
  # If S14 behaved like the null, its max Pi distribution would sit on
  # top of S11c's. AUC = P(max Pi in S14 > max Pi in S11c).
  auc_S14_vs_S11c = mean(outer(d14$max_pi, d11c$max_pi, ">")) +
                    0.5 * mean(outer(d14$max_pi, d11c$max_pi, "==")),
  auc_S1c_vs_S11c = mean(outer(d1c$max_pi, d11c$max_pi, ">")) +
                    0.5 * mean(outer(d1c$max_pi, d11c$max_pi, "=="))
)

# ---- per-criterion resolution ---------------------------------------
crit_summ <- NULL
if (!is.null(ctr_all)) {
  agg <- aggregate(delta ~ scenario + criterion, data = ctr_all,
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                       sd = sd(x, na.rm = TRUE),
                                       pos = mean(x > 0, na.rm = TRUE)))
  crit_summ <- data.frame(scenario = agg$scenario, criterion = agg$criterion,
                          mean_delta = round(agg$delta[, "mean"], 4),
                          sd_delta   = round(agg$delta[, "sd"], 4),
                          pct_favours_true = round(agg$delta[, "pos"] * 100, 1))
  crit_summ <- crit_summ[order(crit_summ$scenario, -abs(crit_summ$mean_delta)), ]
}

# ---- write outputs ---------------------------------------------------
write.csv(ev_all,  file.path(OUT, "ext_per_replicate.csv"), row.names = FALSE)
write.csv(summ,    file.path(OUT, "ext_summary.csv"), row.names = FALSE)
if (!is.null(crit_summ))
  write.csv(crit_summ, file.path(OUT, "ext_criteria_contrast.csv"), row.names = FALSE)

write_json(list(
  meta = list(
    n_rep_requested = N_REP,
    generated = as.character(Sys.time()),
    R_version = R.version.string,
    seed_formula = "seed(scenario, rep) = 20260524 + 17*rep + 1000*match(scenario, SCENARIOS)",
    scenario_order = SCENARIOS,
    total_elapsed_s = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
    per_scenario = meta),
  summary = summ,
  congener_similarity_fit = sim_fit,
  dual_share_fit = share_fit,
  specificity_S14 = spec,
  criteria_contrast = crit_summ
), file.path(OUT, "ext_benchmark_results.json"), auto_unbox = TRUE, pretty = TRUE,
   digits = 6, na = "null")

cat("\n--- S12 congener similarity effect ---\n"); str(sim_fit)
cat("\n--- S13b share effect ---\n"); str(share_fit)
cat("\n--- S14 specificity ---\n"); str(spec)
cat("\nWrote to ", OUT, "\n", sep = "")
cat("Total elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units="mins")), 1),
    " min\n", sep = "")
