#' Criterion C1 — Spatial co-distribution
#'
#' Tests whether human seroprevalence is spatially correlated with animal
#' prevalence of each candidate host, after mutual adjustment.
#' @param data spillover_data object.
#' @return numeric vector of length K, scores on [0,1] (softmax-normalized).
#' @export
criterion_C1_spatial <- function(data) {
  K <- data$n_hosts
  ph <- aggregate(sero ~ site, data = data$human, FUN = mean)
  names(ph)[2] <- "p_human"
  d <- merge(ph, data$site, by = "site")
  prev_cols <- paste0("prev_pcr_h", seq_len(K))
  miss <- setdiff(prev_cols, names(d))
  if (length(miss) > 0)
    return(stats::setNames(rep(0, K), data$host_names))
  formula_str <- paste("p_human ~", paste(prev_cols, collapse = " + "))
  if ("anthrop" %in% names(d)) formula_str <- paste(formula_str, "+ anthrop")
  fit <- tryCatch(stats::lm(as.formula(formula_str), data = d), error = function(e) NULL)
  if (is.null(fit)) return(stats::setNames(rep(0, K), data$host_names))
  cf <- summary(fit)$coefficients
  z <- sapply(prev_cols, function(v) if (v %in% rownames(cf)) cf[v, "t value"] else 0)
  stats::setNames(softmax_norm(z), data$host_names)
}

#' Criterion C2 — Reservoir versus accidental host signatures
#'
#' Combines (i) spatial stability of prevalence and (ii) age-prevalence slope
#' in the animal population to discriminate reservoirs (stable + age-graded)
#' from accidental hosts.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C2_reservoir_signatures <- function(data) {
  K <- data$n_hosts
  cv <- sapply(seq_len(K), function(k) {
    pv <- data$animal$prev_sero[data$animal$host == k]
    sd(pv) / (mean(pv) + 0.01)
  })
  stability <- 1 - rank_norm(cv)
  age_slope <- sapply(seq_len(K), function(k) {
    if (is.null(data$animal_indiv)) return(0)
    d <- data$animal_indiv[data$animal_indiv$host == k, ]
    if (nrow(d) < 20 || length(unique(d$sero)) < 2) return(0)
    f <- tryCatch(stats::glm(sero ~ age, data = d, family = stats::binomial()),
                  error = function(e) NULL)
    if (is.null(f)) 0 else as.numeric(stats::coef(f)["age"])
  })
  age_norm <- rank_norm(age_slope)
  out <- 0.5 * stability + 0.5 * age_norm
  stats::setNames(out, data$host_names)
}

#' Criterion C3 — Individual-level exposure association
#'
#' Multi-level logistic regression: human seropositivity on exposure to each
#' candidate host, with site random effect.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C3_individual <- function(data) {
  K <- data$n_hosts
  d <- data$human
  exp_cols <- paste0("exp_h", seq_len(K))
  covariates <- intersect(c("age", "sex", "prof", "env_declared"), names(d))
  rhs <- paste(c(exp_cols, covariates), collapse = " + ")
  formula_str <- paste("sero ~", rhs, "+ (1 | site)")
  fit <- tryCatch(
    lme4::glmer(stats::as.formula(formula_str), data = d, family = stats::binomial(),
                control = lme4::glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (.glmer_ok(fit, exp_cols)) {
    cf <- summary(fit)$coefficients
    z <- sapply(exp_cols, function(v) if (v %in% rownames(cf)) cf[v, "z value"] else 0)
    return(stats::setNames(softmax_norm(z), data$host_names))
  }
  # Fall back to a fixed-effect site model
  formula_str <- paste("sero ~", rhs, "+ factor(site)")
  fit <- tryCatch(stats::glm(stats::as.formula(formula_str), data = d,
                             family = stats::binomial()),
                  error = function(e) NULL, warning = function(w) NULL)
  if (!.glm_ok(fit, exp_cols))
    return(stats::setNames(rep(NA_real_, K), data$host_names))
  cf <- summary(fit)$coefficients
  z <- sapply(exp_cols, function(v) if (v %in% rownames(cf)) cf[v, "z value"] else 0)
  stats::setNames(softmax_norm(z), data$host_names)
}

#' Criterion C4 — Dose-response (residence time)
#'
#' Tests increasing risk with residence time among individuals declaring
#' exposure to host k.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C4_dose_response <- function(data) {
  K <- data$n_hosts
  d <- data$human
  if (!"residence" %in% names(d)) return(stats::setNames(rep(0, K), data$host_names))
  slopes <- sapply(seq_len(K), function(k) {
    exp_var <- paste0("exp_h", k)
    sub <- d[d[[exp_var]] == 1, ]
    if (nrow(sub) < 20 || length(unique(sub$sero)) < 2) return(0)
    f <- tryCatch(stats::glm(sero ~ residence, data = sub, family = stats::binomial()),
                  error = function(e) NULL)
    if (is.null(f)) 0 else as.numeric(stats::coef(f)["residence"])
  })
  stats::setNames(rank_norm(slopes), data$host_names)
}

#' Criterion C5 — Low cryptic exposure fraction
#'
#' Among seropositives, the fraction explained by declared exposure to host k.
#' Higher score = larger fraction of seropositives are exposed to host k.
#' @param data spillover_data object.
#' @return numeric vector of length K on [0,1].
#' @export
criterion_C5_cryptic_fraction <- function(data) {
  K <- data$n_hosts
  d <- data$human
  pos <- d[d$sero == 1, ]
  if (nrow(pos) == 0) return(stats::setNames(rep(0, K), data$host_names))
  scores <- sapply(seq_len(K), function(k) 1 - mean(pos[[paste0("exp_h", k)]] == 0))
  stats::setNames(rank_norm(scores), data$host_names)
}
