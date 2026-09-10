#' Criterion C6 — Interaction exposure x local prevalence
#'
#' Tests whether the risk associated with exposure to host k is amplified
#' at sites with higher prevalence of host k (signature of active spillover).
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C6_interaction <- function(data) {
  K <- data$n_hosts
  d <- data$human
  for (k in seq_len(K)) {
    prev_col <- paste0("prev_pcr_h", k)
    if (!prev_col %in% names(data$site)) next
    d[[paste0("prev_local_h", k)]] <- data$site[[prev_col]][match(d$site, data$site$site)]
  }
  covariates <- intersect(c("age", "sex", "prof"), names(d))
  z_int <- sapply(seq_len(K), function(k) {
    exp_var  <- paste0("exp_h", k)
    prev_var <- paste0("prev_local_h", k)
    if (!prev_var %in% names(d)) return(NA_real_)
    rhs <- paste0(exp_var, " * ", prev_var)
    if (length(covariates) > 0) rhs <- paste(rhs, "+", paste(covariates, collapse = " + "))
    formula_str <- paste("sero ~", rhs, "+ factor(site)")
    fit <- tryCatch(stats::glm(stats::as.formula(formula_str), data = d,
                               family = stats::binomial()),
                    error = function(e) NULL, warning = function(w) NULL)
    int_term <- paste0(exp_var, ":", prev_var)
    if (!.glm_ok(fit, int_term)) return(NA_real_)
    cf <- summary(fit)$coefficients
    if (int_term %in% rownames(cf)) cf[int_term, "z value"] else 0
  })
  # The scores are normalized across hosts, so a partially estimated vector is
  # not comparable: if any host's interaction fit failed, the criterion is NA.
  if (any(!is.finite(z_int)))
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  stats::setNames(softmax_norm(z_int), data$host_names)
}

#' Criterion C7 — Catalytic FOI coherence
#'
#' Estimates the force of infection in humans (catalytic model on age) and
#' in each candidate host. Scores hosts by raw animal FOI when the ratio
#' rho = lambda_H / lambda_k is biologically plausible (0.001 < rho < 1).
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C7_foi_coherence <- function(data) {
  K <- data$n_hosts
  d_h <- data$human
  if (sum(d_h$sero) < 5) return(stats::setNames(rep(0, K), data$host_names))
  # Estimate lambda_H by simple catalytic via nls
  lam_h <- tryCatch({
    f <- stats::nls(sero ~ 1 - exp(-lambda * age), data = d_h,
                    start = list(lambda = 0.01),
                    control = stats::nls.control(maxiter = 100, warnOnly = TRUE))
    as.numeric(stats::coef(f)["lambda"])
  }, error = function(e) NA)
  if (is.na(lam_h) || lam_h <= 0)
    return(stats::setNames(rep(0, K), data$host_names))

  scores <- sapply(seq_len(K), function(k) {
    if (is.null(data$animal_indiv)) return(0)
    d_a <- data$animal_indiv[data$animal_indiv$host == k, ]
    if (sum(d_a$sero) < 5 || length(unique(d_a$sero)) < 2) return(0)
    lam_k <- tryCatch({
      f <- stats::nls(sero ~ 1 - exp(-lambda * age), data = d_a,
                      start = list(lambda = 0.5),
                      control = stats::nls.control(maxiter = 100, warnOnly = TRUE))
      as.numeric(stats::coef(f)["lambda"])
    }, error = function(e) NA)
    if (is.na(lam_k) || lam_k <= 0) return(0)
    rho <- lam_h / lam_k
    if (rho > 0.001 && rho < 1) lam_k else lam_k * 0.5
  })
  stats::setNames(rank_norm(scores), data$host_names)
}

#' Criterion C8 — Temporal precedence
#'
#' Cross-correlation between animal prevalence time series and human incidence
#' time series, at non-negative lags (animal precedes human).
#' Returns 0 for all hosts when temporal data are missing.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1], or NA if not applicable.
#' @export
criterion_C8_temporal_precedence <- function(data) {
  K <- data$n_hosts
  if (is.null(data$animal_ts) || is.null(data$human_ts))
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  h_yearly <- aggregate(incidence ~ year, data = data$human_ts, FUN = mean)
  scores <- sapply(seq_len(K), function(k) {
    a_k <- data$animal_ts[data$animal_ts$host == k, ]
    if (nrow(a_k) < 4) return(0)
    a_k$year_floor <- floor(a_k$year)
    a_yearly <- aggregate(prev_pcr_true ~ year_floor, data = a_k, FUN = mean)
    names(a_yearly)[1] <- "year"
    m <- merge(h_yearly, a_yearly, by = "year")
    if (nrow(m) < 5) return(0)
    best <- -1
    for (lag in 0:3) {
      x <- m$prev_pcr_true
      y <- m$incidence
      if (lag > 0) {
        x <- utils::head(x, -lag)
        y <- utils::tail(y, length(x))
      }
      if (length(x) < 4 || sd(x) == 0 || sd(y) == 0) next
      r <- stats::cor(x, y)
      if (is.finite(r) && r > best) best <- r
    }
    max(best, 0)
  })
  stats::setNames(rank_norm(scores), data$host_names)
}

#' Criterion C9 — Mediation by animal prevalence
#'
#' Difference-of-coefficients mediation: how much of the anthropization
#' effect on human seropositivity transits through the prevalence of host k.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C9_mediation <- function(data) {
  K <- data$n_hosts
  d <- data$human
  if (!"anthrop" %in% names(d))
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  for (k in seq_len(K)) {
    prev_col <- paste0("prev_pcr_h", k)
    if (prev_col %in% names(data$site)) {
      d[[paste0("prev_local_h", k)]] <- data$site[[prev_col]][match(d$site, data$site$site)]
    }
  }
  covariates <- intersect(c("age", "sex", "prof"), names(d))
  rhs_base <- if (length(covariates) > 0) paste("anthrop +", paste(covariates, collapse = " + "))
              else "anthrop"
  fit_total <- tryCatch(stats::glm(stats::as.formula(paste("sero ~", rhs_base)),
                                   data = d, family = stats::binomial()),
                        error = function(e) NULL)
  if (is.null(fit_total)) return(stats::setNames(rep(0, K), data$host_names))
  b_total <- stats::coef(fit_total)["anthrop"]
  if (is.na(b_total) || b_total <= 0)
    return(stats::setNames(rep(0, K), data$host_names))
  scores <- sapply(seq_len(K), function(k) {
    prev_var <- paste0("prev_local_h", k)
    if (!prev_var %in% names(d)) return(0)
    rhs <- paste(rhs_base, "+", prev_var)
    fit_dir <- tryCatch(stats::glm(stats::as.formula(paste("sero ~", rhs)),
                                   data = d, family = stats::binomial()),
                        error = function(e) NULL)
    if (is.null(fit_dir)) return(0)
    b_dir <- stats::coef(fit_dir)["anthrop"]
    max((b_total - b_dir) / b_total, 0)
  })
  stats::setNames(rank_norm(scores), data$host_names)
}

#' Criterion C10 — Mechanistic Akaike weights
#'
#' Fits one mechanistic model per candidate source-k and computes Akaike
#' weights across models.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1] summing to 1.
#' @export
criterion_C10_mechanistic <- function(data) {
  K <- data$n_hosts
  d <- data$human
  for (k in seq_len(K)) {
    prev_col <- paste0("prev_pcr_h", k)
    if (prev_col %in% names(data$site)) {
      d[[paste0("prev_local_h", k)]] <- data$site[[prev_col]][match(d$site, data$site$site)]
    }
  }
  covariates <- intersect(c("age", "sex", "prof", "env_declared", "residence"), names(d))
  aics <- sapply(seq_len(K), function(k) {
    exp_var  <- paste0("exp_h", k)
    prev_var <- paste0("prev_local_h", k)
    if (!prev_var %in% names(d)) return(Inf)
    rhs <- paste0("I(", exp_var, " * ", prev_var, ")")
    if (length(covariates) > 0) rhs <- paste(rhs, "+", paste(covariates, collapse = " + "))
    fit <- tryCatch(stats::glm(stats::as.formula(paste("sero ~", rhs)),
                               data = d, family = stats::binomial()),
                    error = function(e) NULL)
    if (is.null(fit)) Inf else stats::AIC(fit)
  })
  if (all(is.infinite(aics))) return(stats::setNames(rep(0, K), data$host_names))
  delta <- aics - min(aics, na.rm = TRUE)
  w <- exp(-0.5 * delta)
  w <- w / sum(w)
  stats::setNames(w, data$host_names)
}

#' Criterion C11 — Negative control specificity
#'
#' Returns NA when no negative control data are provided. When provided
#' (\code{data$controls}), checks that signal is absent for the control
#' pathogen / exposure / host.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1] or NA.
#' @export
criterion_C11_specificity <- function(data) {
  K <- data$n_hosts
  if (is.null(data$controls))
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  # If controls provided, simple implementation: score = 1 if control C3
  # for that host is NOT significant, otherwise 0.
  ctrl <- data$controls
  if (is.null(ctrl$control_human)) return(stats::setNames(rep(NA_real_, K), data$host_names))
  d_ctrl <- ctrl$control_human
  scores <- sapply(seq_len(K), function(k) {
    if (length(unique(d_ctrl$sero)) < 2) return(1)  # no variation = trivially no signal
    f <- tryCatch(stats::glm(stats::as.formula(paste0("sero ~ exp_h", k)),
                             data = d_ctrl, family = stats::binomial()),
                  error = function(e) NULL)
    if (is.null(f)) return(0.5)
    p <- summary(f)$coefficients[paste0("exp_h", k), "Pr(>|z|)"]
    if (is.na(p)) 0.5 else if (p > 0.05) 1 else 0
  })
  stats::setNames(scores, data$host_names)
}
