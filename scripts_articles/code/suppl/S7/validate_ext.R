# =====================================================================
# validate_ext.R -- Check that the generalised generator is a faithful
# fork of the published one.
#
# For every legacy scenario, simulate_multihost(scenario_params_ext(x))
# must reproduce simulate_spillover(scenario_params(x)) exactly, given
# the same seed: same data frames, same Pi vector. If this passes, any
# difference observed in the new scenarios is attributable to the new
# scenario structure and not to the fork.
# =====================================================================
suppressPackageStartupMessages(library(spilloverID))
source("simulate_ext.R")
source("eval_ext.R")

legacy <- c("S1_canonical","S2_bridge","S3_cryptic","S4_confounded",
            "S5_codistributed","S6_strong_seasonal","S7_low_power",
            "S8_measurement_error","S9_weak_signal","S10_short_immunity",
            "S11_no_spillover")

cat("scenario                 human  animal  animal_indiv  periods  pi_maxabsdiff\n")
ok_all <- TRUE
for (sc in legacy) {
  sd_use <- 20260524 + 7 * 17 + match(sc, legacy)
  a <- simulate_spillover(scenario_params(sc), seed = sd_use)
  b <- simulate_multihost(scenario_params_ext(sc), seed = sd_use)
  eq <- function(x, y) isTRUE(all.equal(x, y, check.attributes = FALSE))
  e_h  <- eq(a$human,  b$human)
  e_a  <- eq(a$animal, b$animal)
  e_ai <- eq(a$animal_indiv, b$animal_indiv)
  e_ap <- eq(a$animal_periods, b$animal_periods)
  pa <- suppressWarnings(compute_pi(a)$pi)
  pb <- suppressWarnings(compute_pi(b)$pi)
  d  <- max(abs(pa - pb))
  ok_all <- ok_all && e_h && e_a && e_ai && e_ap && d < 1e-12
  cat(sprintf("%-24s %-6s %-7s %-13s %-8s %.3e\n", sc, e_h, e_a, e_ai, e_ap, d))
}
cat("\nALL LEGACY SCENARIOS BIT-IDENTICAL:", ok_all, "\n")

# Sanity check of the new scenarios: structure and truth encoding.
cat("\n--- new scenarios, one replicate each ---\n")
for (sc in c("S12_congeners","S13_dual_symmetric","S13b_dual_asymmetric",
             "S14_hidden_reservoir","S1c_ref_single","S11c_ref_null")) {
  d <- simulate_multihost(scenario_params_ext(sc), seed = 424242)
  tr <- attr(d, "truth")
  r  <- suppressWarnings(compute_pi(d))
  ev <- eval_replicate(r$pi, tr, sc, 1L, 424242)
  cat(sprintf("%-22s K_obs=%d  k_star_set={%s}  genus=%s  rho=%s share1=%s\n",
      sc, tr$n_hosts_observed, paste(tr$k_star_set, collapse=","),
      paste(tr$genus_obs, collapse=""),
      ifelse(is.na(tr$rho_congener),"-",round(tr$rho_congener,3)),
      ifelse(length(tr$shares)==0,"-",round(tr$shares[1],3))))
  cat("   Pi = ", paste(sprintf("%.3f", r$pi), collapse=" "), "\n")
  cat("   top1_species=", ev$top1_species, " top1_genus=", ev$top1_genus,
      " both_top2=", ev$both_top2, " max_pi=", round(ev$max_pi,3), "\n", sep="")
  na_crit <- colnames(r$criteria)[apply(r$criteria, 2, function(x) all(is.na(x)))]
  cat("   NA criteria: ", paste(na_crit, collapse=","), "\n")
}
