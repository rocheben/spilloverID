#' Rank-based normalization
#'
#' Maps a vector to evenly spaced values in [0,1] based on ranks.
#' @param x numeric vector.
#' @return numeric vector of same length on [0,1].
#' @export
rank_norm <- function(x) {
  x <- as.numeric(x)
  K <- length(x)
  if (K <= 1) return(rep(1, K))
  (rank(x) - 1) / (K - 1)
}

#' Softmax-style normalization
#'
#' Standardize the input, apply softmax with given temperature, then rescale
#' so the maximum value is 1. Reflects the magnitude of evidence rather than
#' only ranks.
#' @param x numeric vector.
#' @param temperature softmax temperature (default 1).
#' @return numeric vector on [0,1] with max = 1.
#' @export
softmax_norm <- function(x, temperature = 1) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  if (all(x == 0)) return(rep(0, length(x)))
  x_std <- (x - mean(x)) / (sd(x) + 1e-6)
  e <- exp(x_std / temperature)
  p <- e / sum(e)
  p / max(p)
}

#' Default weights for the composite Pi score
#'
#' Returns the preset weights for criteria C1 to C14.
#' @return named numeric vector summing to 1.
#' @export
default_weights <- function() {
  # Original C1-C11 weights from Roche et al. (initial framework)
  # contracted by factor 0.55; C12-C14 (maintenance signatures) carry 0.45.
  # This allocation reflects the empirical finding that the maintenance vs
  # amplification distinction is the dominant axis of inferential error
  # in retrospective outbreak analysis (see Methods).
  w <- c(C1 = 0.04, C2 = 0.04, C3 = 0.06, C4 = 0.04, C5 = 0.03,
         C6 = 0.09, C7 = 0.04, C8 = 0.03, C9 = 0.06, C10 = 0.07,
         C11 = 0.05, C12 = 0.16, C13 = 0.14, C14 = 0.15)
  w / sum(w)
}


# --- Convergence diagnostics -------------------------------------------------
# C3 and C6 rely on mixed-effects and interaction fits that can fail to converge
# on sparse reconstructions (degenerate Hessian in glmer, complete or
# quasi-complete separation in glm). A non-converged fit, or a separated term of
# interest, yields a coefficient that is numerically unstable between runs, so
# the criterion returns NA and compute_pi() redistributes its weight rather than
# propagating a drifting value. Separation is diagnosed on the terms the
# criterion actually reads, not on the nuisance site effects: a site in which no
# individual seroconverted drives its own indicator to a boundary without making
# the exposure coefficient unidentifiable. These helpers are internal.

.terms_estimable <- function(cf, terms, max_abs_coef = 30, max_se = 30) {
  present <- intersect(terms, rownames(cf))
  if (length(present) == 0) return(TRUE)
  est <- cf[present, 1]
  se  <- cf[present, 2]
  if (any(!is.finite(est)) || any(!is.finite(se))) return(FALSE)
  if (any(abs(est) > max_abs_coef) || any(se > max_se)) return(FALSE)
  TRUE
}

.glmer_ok <- function(fit, terms) {
  if (is.null(fit)) return(FALSE)
  msgs <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  if (length(msgs) > 0) return(FALSE)
  h <- tryCatch(fit@optinfo$derivs$Hessian, error = function(e) NULL)
  if (!is.null(h)) {
    ev <- tryCatch(eigen(h, only.values = TRUE)$values, error = function(e) NA_real_)
    if (any(!is.finite(ev)) || min(Re(ev)) <= 0) return(FALSE)
  }
  cf <- tryCatch(stats::coef(summary(fit)), error = function(e) NULL)
  if (is.null(cf)) return(FALSE)
  .terms_estimable(cf, terms)
}

.glm_ok <- function(fit, terms) {
  if (is.null(fit)) return(FALSE)
  if (!isTRUE(fit$converged)) return(FALSE)
  cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(cf)) return(FALSE)
  .terms_estimable(cf, terms)
}
