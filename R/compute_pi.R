#' Compute the composite plausibility score Pi for each candidate host
#'
#' Calls each of the 14 criterion functions on the data, normalizes scores
#' to [0,1], and returns a weighted sum. Criteria returning NA (e.g. C8
#' without temporal data, C11 without controls, C12 without paired
#' epidemic/inter-epidemic samples) have their weight redistributed
#' proportionally over the active criteria.
#'
#' @param data spillover_data object.
#' @param weights named numeric vector with names C1..C14. Defaults to
#'   \code{\link{default_weights}}.
#' @param bootstrap integer; number of site-level bootstrap replicates
#'   for confidence intervals. Default 0 (no bootstrap).
#' @param verbose logical; print progress.
#'
#' @return A list with:
#'   \item{pi}{named numeric vector of Pi scores per host.}
#'   \item{criteria}{matrix of per-criterion scores (rows = hosts, cols = criteria).}
#'   \item{weights_used}{actual weights applied (after redistribution).}
#'   \item{rank}{integer ranks of hosts (1 = highest Pi).}
#'   \item{bootstrap}{matrix of bootstrap Pi if requested, NULL otherwise.}
#' @export
compute_pi <- function(data, weights = NULL, bootstrap = 0, verbose = FALSE) {
  if (is.null(weights)) weights <- default_weights()
  weights <- weights[names(weights) %in% paste0("C", 1:14)]
  weights <- weights / sum(weights)

  crit_funs <- list(
    C1 = criterion_C1_spatial, C2 = criterion_C2_reservoir_signatures,
    C3 = criterion_C3_individual, C4 = criterion_C4_dose_response,
    C5 = criterion_C5_cryptic_fraction, C6 = criterion_C6_interaction,
    C7 = criterion_C7_foi_coherence, C8 = criterion_C8_temporal_precedence,
    C9 = criterion_C9_mediation, C10 = criterion_C10_mechanistic,
    C11 = criterion_C11_specificity,
    C12 = criterion_C12_inter_epidemic,
    C13 = criterion_C13_endemic_age,
    C14 = criterion_C14_spatial_independence
  )
  K <- data$n_hosts
  crit_mat <- matrix(NA_real_, nrow = K, ncol = length(crit_funs),
                     dimnames = list(data$host_names, names(crit_funs)))
  for (cn in names(crit_funs)) {
    if (verbose) message("  Computing ", cn, " ...")
    res <- tryCatch(crit_funs[[cn]](data), error = function(e) rep(NA_real_, K))
    crit_mat[, cn] <- res
  }

  na_crit <- apply(crit_mat, 2, function(x) all(is.na(x)))
  w_active <- weights
  w_active[names(w_active) %in% names(na_crit[na_crit])] <- 0
  if (sum(w_active) > 0) w_active <- w_active / sum(w_active)

  crit_safe <- crit_mat
  crit_safe[is.na(crit_safe)] <- 0
  pi_vec <- as.numeric(crit_safe %*% w_active)
  names(pi_vec) <- data$host_names

  out <- list(
    pi = pi_vec,
    criteria = crit_mat,
    weights_used = w_active,
    rank = rank(-pi_vec, ties.method = "min"),
    bootstrap = NULL
  )

  if (bootstrap > 0) {
    if (verbose) message("Bootstrap (", bootstrap, " replicates)...")
    sites <- unique(data$human$site)
    boot_mat <- matrix(NA_real_, nrow = bootstrap, ncol = K)
    colnames(boot_mat) <- data$host_names
    for (b in seq_len(bootstrap)) {
      sampled <- sample(sites, replace = TRUE)
      d_b <- .resample_data(data, sampled)
      pi_b <- tryCatch(compute_pi(d_b, weights = weights, bootstrap = 0)$pi,
                       error = function(e) rep(NA_real_, K))
      boot_mat[b, ] <- pi_b
    }
    out$bootstrap <- boot_mat
  }
  class(out) <- "spillover_pi"
  out
}

#' @keywords internal
.resample_data <- function(data, sampled_sites) {
  new_human <- list(); new_animal <- list(); new_animal_indiv <- list(); new_site <- list()
  for (i in seq_along(sampled_sites)) {
    s <- sampled_sites[i]; new_id <- i
    h <- data$human[data$human$site == s, ]; h$site <- new_id; new_human[[i]] <- h
    a <- data$animal[data$animal$site == s, ]; a$site <- new_id; new_animal[[i]] <- a
    if (!is.null(data$animal_indiv)) {
      ai <- data$animal_indiv[data$animal_indiv$site == s, ]; ai$site <- new_id
      new_animal_indiv[[i]] <- ai
    }
    sr <- data$site[data$site$site == s, ]; sr$site <- new_id; new_site[[i]] <- sr
  }
  data2 <- data
  data2$human <- do.call(rbind, new_human)
  data2$animal <- do.call(rbind, new_animal)
  if (!is.null(data$animal_indiv)) data2$animal_indiv <- do.call(rbind, new_animal_indiv)
  data2$site <- do.call(rbind, new_site)
  data2
}

#' @export
print.spillover_pi <- function(x, ...) {
  cat("Composite plausibility scores Pi:\n")
  df <- data.frame(host = names(x$pi), Pi = round(x$pi, 3), rank = x$rank)
  df <- df[order(df$rank), ]
  print(df, row.names = FALSE)
  cat("\nCriterion scores:\n")
  print(round(x$criteria, 3))
  cat("\nWeights used (after NA redistribution):\n")
  print(round(x$weights_used, 3))
  if (!is.null(x$bootstrap)) {
    cat("\nBootstrap CI (95%):\n")
    ci <- apply(x$bootstrap, 2, stats::quantile,
                probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    print(round(t(ci), 3))
  }
  invisible(x)
}
