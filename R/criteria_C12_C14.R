#' Criterion C12 — Inter-epidemic persistence
#'
#' Compares seroprevalence in each candidate host between epidemic and
#' inter-epidemic periods. A genuine maintenance reservoir sustains
#' circulation in the absence of human cases (ratio close to 1); an
#' amplifying or bridge host shows transient pulses tied to outbreaks
#' (ratio close to 0).
#'
#' Requires \code{data$animal_periods}, a data frame with columns
#' \code{site}, \code{host}, \code{period} (factor with at least levels
#' "epidemic" and "interepidemic"), \code{n_sampled}, \code{n_sero}.
#'
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1], or NA if data unavailable.
#' @export
criterion_C12_inter_epidemic <- function(data) {
  K <- data$n_hosts
  if (is.null(data$animal_periods)) {
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  }
  ap <- data$animal_periods
  required <- c("host", "period", "n_sampled", "n_sero")
  if (!all(required %in% names(ap))) {
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  }
  ratios <- sapply(seq_len(K), function(k) {
    sub <- ap[ap$host == k, ]
    if (nrow(sub) == 0) return(0)
    epi <- sub[sub$period == "epidemic", ]
    inter <- sub[sub$period == "interepidemic", ]
    if (nrow(epi) == 0 || nrow(inter) == 0) return(0)
    p_epi <- sum(epi$n_sero) / max(sum(epi$n_sampled), 1)
    p_inter <- sum(inter$n_sero) / max(sum(inter$n_sampled), 1)
    if (p_epi <= 0) return(0)
    min(p_inter / p_epi, 1)
  })
  stats::setNames(softmax_norm(ratios), data$host_names)
}

#' Criterion C13 — Endemic age-prevalence signature
#'
#' Distinguishes endemic (cumulative-with-age, catalytic) from epidemic
#' (age-independent, cohort) seropositivity patterns. For each host, fits
#' a catalytic model and a constant-prevalence model, then computes the
#' Akaike weight for the catalytic (endemic) model.
#'
#' Requires \code{data$animal_indiv} with individual-level \code{age}, \code{sero},
#' and \code{host} columns.
#'
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C13_endemic_age <- function(data) {
  K <- data$n_hosts
  if (is.null(data$animal_indiv)) {
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  }
  w_endemic <- sapply(seq_len(K), function(k) {
    d_k <- data$animal_indiv[data$animal_indiv$host == k, ]
    if (nrow(d_k) < 30 || length(unique(d_k$sero)) < 2) return(0.5)
    # Endemic catalytic: logit(P) = b0 + b1 * age  (b1 positive for endemic accumulation)
    fit_endemic <- tryCatch(
      stats::glm(sero ~ age, data = d_k, family = stats::binomial()),
      error = function(e) NULL
    )
    # Cohort/epidemic: constant probability across ages (intercept only)
    fit_constant <- tryCatch(
      stats::glm(sero ~ 1, data = d_k, family = stats::binomial()),
      error = function(e) NULL
    )
    if (is.null(fit_endemic) || is.null(fit_constant)) return(0.5)
    # Penalise negative slope (would indicate antibody waning, not accumulation)
    age_coef <- as.numeric(stats::coef(fit_endemic)["age"])
    if (is.na(age_coef) || age_coef <= 0) return(0.1)
    aic_end <- stats::AIC(fit_endemic)
    aic_con <- stats::AIC(fit_constant)
    delta_end <- aic_end - min(aic_end, aic_con)
    delta_con <- aic_con - min(aic_end, aic_con)
    w <- exp(-0.5 * delta_end) / (exp(-0.5 * delta_end) + exp(-0.5 * delta_con))
    w
  })
  stats::setNames(softmax_norm(w_endemic), data$host_names)
}

#' Criterion C14 — Spatial independence from co-circulating hosts
#'
#' Tests whether each candidate host's seroprevalence is preserved in sites
#' where the OTHER candidate hosts (potential amplifiers) are at lowest
#' prevalence. A genuine maintenance reservoir maintains its own positivity
#' independently of these other hosts; an amplifier dependent on a reservoir
#' for repeated introduction will collapse when its source is absent.
#'
#' For each host k, defines "low-other sites" as those where the maximum
#' prevalence among all OTHER hosts is below its median across sites. Computes
#' the ratio of mean prevalence in low-other sites to overall mean prevalence
#' for host k.
#'
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C14_spatial_independence <- function(data) {
  K <- data$n_hosts
  if (K < 2) return(stats::setNames(rep(NA_real_, K), data$host_names))
  # Build prevalence matrix sites x hosts
  S <- nrow(data$site)
  P <- matrix(NA_real_, nrow = S, ncol = K)
  for (k in seq_len(K)) {
    col_k <- paste0("prev_sero_h", k)
    if (col_k %in% names(data$site)) P[, k] <- data$site[[col_k]]
  }
  if (any(is.na(P))) P[is.na(P)] <- 0
  scores <- sapply(seq_len(K), function(k) {
    others_idx <- setdiff(seq_len(K), k)
    if (length(others_idx) == 0) return(0)
    max_others <- apply(P[, others_idx, drop = FALSE], 1, max, na.rm = TRUE)
    med_max_others <- stats::median(max_others, na.rm = TRUE)
    low_other_sites <- max_others <= med_max_others
    if (sum(low_other_sites, na.rm = TRUE) < 2) return(0)
    p_k_lowother <- mean(P[low_other_sites, k], na.rm = TRUE)
    p_k_all <- mean(P[, k], na.rm = TRUE)
    # Apply prevalence floor to avoid trivial high score for hosts with
    # near-zero prevalence everywhere
    if (p_k_all <= 0.01) return(0)
    # Multiply by mean prevalence to favour hosts that ARE present
    # (so a host with low overall prevalence does not get high score just
    # because it is rare everywhere)
    stability <- min(p_k_lowother / p_k_all, 1)
    stability * sqrt(p_k_all)  # weight by sqrt of overall prevalence
  })
  stats::setNames(softmax_norm(scores), data$host_names)
}
