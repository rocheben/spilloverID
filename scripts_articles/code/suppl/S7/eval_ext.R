# =====================================================================
# eval_ext.R -- Evaluation metrics for scenarios whose ground truth is
# not a scalar k_star.
#
# run_benchmark.R assumes `k_star <- attr(d, "truth")$k_star` is a single
# integer, and scores top-1 / top-2 against it. That assumption fails for
# three of the new scenarios:
#   * S12  the truth is a species WITHIN a genus, so "correct" has two
#          nested meanings (exact species vs correct genus);
#   * S13  the truth is a SET of two hosts, so top-1 is not even the
#          right question;
#   * S14  the truth is not in the candidate set at all, so every rank
#          is meaningless and only the calibration of max(Pi) matters.
# This module returns one wide row per replicate covering all three
# regimes, with NA in the fields that do not apply. Legacy scenarios pass
# through it unchanged: for a scalar truth, top1_species and top1_genus
# both reduce to the published top-1 definition.
# =====================================================================

eval_replicate <- function(pi_vec, truth, scenario, rep, seed) {
  K <- length(pi_vec)
  rk <- rank(-pi_vec, ties.method = "min")
  ks <- truth$k_star_set
  genus <- truth$genus_obs
  ord <- order(-pi_vec)
  top1 <- ord[1]; top2_set <- ord[seq_len(min(2, K))]

  row <- data.frame(
    scenario = scenario, rep = rep, seed = seed,
    n_hosts = K,
    n_true_reservoirs = length(ks),
    reservoir_hidden = isTRUE(truth$reservoir_hidden),
    rho_congener = if (is.null(truth$rho_congener)) NA_real_ else truth$rho_congener,
    share_1 = if (length(truth$shares)) truth$shares[1] else NA_real_,
    top1_host = top1,
    max_pi = max(pi_vec, na.rm = TRUE),
    min_pi = min(pi_vec, na.rm = TRUE),
    pi_range = diff(range(pi_vec, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  for (k in seq_len(K)) row[[paste0("pi_h", k)]] <- pi_vec[k]

  # ---------------- no true reservoir among the candidates -----------
  # (negative control S11, and incomplete inventory S14)
  if (length(ks) == 0) {
    row$top1_species <- NA_integer_
    row$top1_genus   <- NA_integer_
    row$rank_true    <- NA_integer_
    row$pi_true      <- NA_real_
    row$pi_max_other <- max(pi_vec, na.rm = TRUE)
    row$margin       <- NA_real_
    row$both_top2    <- NA_integer_
    row$any_top1     <- NA_integer_
    row$mean_rank_true <- NA_real_
    row$pi_gap_within_true <- NA_real_
    row$margin_set   <- NA_real_
    row$exceeds_055  <- as.integer(max(pi_vec, na.rm = TRUE) > 0.55)
    row$pi_congener  <- NA_real_
    row$rank_congener <- NA_integer_
    return(row)
  }

  # ---------------- single true reservoir ----------------------------
  if (length(ks) == 1) {
    k1 <- ks[1]
    others <- setdiff(seq_len(K), k1)
    row$rank_true    <- rk[k1]
    row$pi_true      <- pi_vec[k1]
    row$pi_max_other <- if (length(others)) max(pi_vec[others]) else NA_real_
    row$margin       <- row$pi_true - row$pi_max_other
    row$top1_species <- as.integer(rk[k1] == 1)
    row$top1_genus   <- as.integer(genus[top1] == genus[k1])
    # congener = any other observed host of the same genus
    cong <- setdiff(which(genus == genus[k1]), k1)
    if (length(cong) == 1) {
      row$pi_congener   <- pi_vec[cong]
      row$rank_congener <- rk[cong]
    } else {
      row$pi_congener   <- NA_real_
      row$rank_congener <- NA_integer_
    }
    row$both_top2    <- NA_integer_
    row$any_top1     <- as.integer(rk[k1] == 1)
    row$mean_rank_true <- as.numeric(rk[k1])
    row$pi_gap_within_true <- NA_real_
    row$margin_set   <- row$margin
    row$exceeds_055  <- as.integer(max(pi_vec, na.rm = TRUE) > 0.55)
    return(row)
  }

  # ---------------- two or more true reservoirs ----------------------
  others <- setdiff(seq_len(K), ks)
  pi_true_set <- pi_vec[ks]
  row$rank_true      <- NA_integer_          # undefined: no unique truth
  row$pi_true        <- max(pi_true_set)     # best-supported true reservoir
  row$pi_max_other   <- if (length(others)) max(pi_vec[others]) else NA_real_
  row$margin         <- row$pi_true - row$pi_max_other
  # recall@|ks|: are ALL true reservoirs in the top |ks|?
  row$both_top2      <- as.integer(all(ks %in% top2_set))
  row$any_top1       <- as.integer(top1 %in% ks)
  row$mean_rank_true <- mean(rk[ks])
  # within-truth spread: large gap => the score arbitrarily designates one
  row$pi_gap_within_true <- max(pi_true_set) - min(pi_true_set)
  # set margin: does the WEAKEST true reservoir still beat every decoy?
  row$margin_set     <- min(pi_true_set) - row$pi_max_other
  row$top1_species   <- as.integer(top1 %in% ks)
  row$top1_genus     <- as.integer(genus[top1] %in% genus[ks])
  row$pi_congener    <- NA_real_
  row$rank_congener  <- NA_integer_
  row$exceeds_055    <- as.integer(max(pi_vec, na.rm = TRUE) > 0.55)
  row
}

# Per-criterion contrast between two hosts, averaged later over replicates.
# Answers: which of C1-C14 actually carries species-level resolution?
criteria_contrast <- function(crit_mat, ka, kb, scenario, rep) {
  cn <- colnames(crit_mat)
  data.frame(scenario = scenario, rep = rep, criterion = cn,
             score_a = as.numeric(crit_mat[ka, ]),
             score_b = as.numeric(crit_mat[kb, ]),
             delta   = as.numeric(crit_mat[ka, ]) - as.numeric(crit_mat[kb, ]),
             stringsAsFactors = FALSE)
}
