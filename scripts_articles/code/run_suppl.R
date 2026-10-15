#' Supplementary analyses for Roche et al. Nature Microbiology
#'
#' Produces:
#'   suppl_sensitivity.json  - weight perturbation sensitivity for 3 cases
#'   suppl_loco.json         - leave-one-criterion-out for 3 cases

suppressPackageStartupMessages({
  library(spilloverID)
  library(jsonlite)
})

set.seed(20260524)
B_SENS <- 200  # number of weight perturbations

# Reload the three case datasets directly from the existing RDS file
cases <- readRDS("../data/case_results.rds")

# ---------------------------------------------------------------------------
# 1. Weight perturbation sensitivity
#    For each case, draw B random Dirichlet perturbations of the weights,
#    compute Pi, and record (a) how often the true host stays at rank 1,
#    (b) the distribution of Pi_true.
# ---------------------------------------------------------------------------
default_w <- default_weights()

sens_results <- list()
for (case_key in c("sn", "nv", "rv")) {
  cc <- cases[[case_key]]
  case_name <- cc$case
  true_host <- cc$true_host
  d <- cc$data
  hosts <- names(cc$result$pi)
  truth_idx <- which(hosts == true_host)

  # Baseline
  base_pi <- cc$result$pi[truth_idx]
  base_rank <- cc$result$rank[truth_idx]

  rank_history <- integer(B_SENS)
  pi_history <- numeric(B_SENS)
  for (b in seq_len(B_SENS)) {
    # Dirichlet perturbation: alpha = default_w * 30 keeps near baseline but adds noise
    alphas <- default_w * 30
    g <- rgamma(length(alphas), shape = alphas, rate = 1)
    w_new <- g / sum(g)
    names(w_new) <- names(default_w)
    res_b <- suppressWarnings(compute_pi(d, weights = w_new, bootstrap = 0))
    rank_history[b] <- res_b$rank[truth_idx]
    pi_history[b] <- res_b$pi[truth_idx]
  }

  sens_results[[case_key]] <- list(
    case = case_name,
    true_host = true_host,
    baseline_pi = base_pi,
    baseline_rank = base_rank,
    perturbation_n = B_SENS,
    pct_rank1 = mean(rank_history == 1) * 100,
    pct_rank_top2 = mean(rank_history <= 2) * 100,
    pi_median = median(pi_history),
    pi_2.5 = quantile(pi_history, 0.025),
    pi_97.5 = quantile(pi_history, 0.975)
  )
  cat(case_name, ": rank-1 stable in",
      round(mean(rank_history == 1) * 100), "% of weight perturbations\n")
}
write_json(sens_results, "../data/suppl_sensitivity.json",
           auto_unbox = TRUE, pretty = TRUE)


# ---------------------------------------------------------------------------
# 2. Leave-one-criterion-out (LOCO)
#    For each case, drop each criterion (set its weight to 0 and redistribute),
#    recompute Pi, and record the resulting rank of the true host and its Pi.
# ---------------------------------------------------------------------------
crit_names <- paste0("C", 1:14)
loco_results <- list()
for (case_key in c("sn", "nv", "rv")) {
  cc <- cases[[case_key]]
  d <- cc$data
  hosts <- names(cc$result$pi)
  truth_idx <- which(hosts == cc$true_host)

  rows <- list()
  for (cn in crit_names) {
    w_new <- default_w
    w_new[cn] <- 0
    w_new <- w_new / sum(w_new)
    res_b <- suppressWarnings(compute_pi(d, weights = w_new, bootstrap = 0))
    rows[[cn]] <- data.frame(
      criterion = cn,
      pi_true = res_b$pi[truth_idx],
      rank_true = res_b$rank[truth_idx],
      delta_pi = res_b$pi[truth_idx] - cc$result$pi[truth_idx]
    )
  }
  loco_results[[case_key]] <- list(
    case = cc$case, true_host = cc$true_host,
    baseline_pi = cc$result$pi[truth_idx],
    baseline_rank = cc$result$rank[truth_idx],
    rows = do.call(rbind, rows)
  )
}
write_json(loco_results, "../data/suppl_loco.json",
           auto_unbox = TRUE, pretty = TRUE)
cat("\nWrote suppl_sensitivity.json and suppl_loco.json\n")
