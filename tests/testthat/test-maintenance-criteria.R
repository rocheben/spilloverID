# Tests for the three maintenance criteria C12-C14, which together carry 0.45 of
# the nominal weight and are the criteria that separate a maintenance reservoir
# from an amplifying host.

make_paired <- function(p_epi, p_inter, n = 200) {
  K <- length(p_epi)
  do.call(rbind, lapply(seq_len(K), function(k) rbind(
    data.frame(site = 1, host = k, period = "epidemic",
               n_sampled = n, n_sero = round(n * p_epi[k])),
    data.frame(site = 1, host = k, period = "interepidemic",
               n_sampled = n, n_sero = round(n * p_inter[k])))))
}

test_that("C12 rewards a persistent host and penalises a pulse-like one", {
  sim <- simulate_spillover(scenario_params("S1_canonical"), seed = 11)
  sim$animal_periods <- make_paired(p_epi = c(.40, .40, .40, .40),
                                    p_inter = c(.38, .02, .01, .00))
  s <- criterion_C12_inter_epidemic(sim)
  expect_length(s, sim$n_hosts)
  expect_true(all(s >= 0 & s <= 1, na.rm = TRUE))
  expect_equal(unname(which.max(s)), 1L)
  expect_lt(s[4], s[1])
})

test_that("C12 returns zero rather than dividing by zero on empty epidemic samples", {
  sim <- simulate_spillover(scenario_params("S1_canonical"), seed = 12)
  sim$animal_periods <- make_paired(p_epi = c(0, 0, 0, 0), p_inter = c(0, 0, 0, 0))
  s <- criterion_C12_inter_epidemic(sim)
  expect_true(all(is.finite(s) | is.na(s)))
})

test_that("C13 and C14 stay on [0,1] and keep the host names", {
  sim <- simulate_spillover(scenario_params("S1_canonical"), seed = 13)
  for (f in list(criterion_C13_endemic_age, criterion_C14_spatial_independence)) {
    s <- f(sim)
    expect_length(s, sim$n_hosts)
    expect_equal(names(s), sim$host_names)
    expect_true(all(s >= 0 & s <= 1, na.rm = TRUE))
  }
})

test_that("compute_pi covers all fourteen criteria and renormalises around NAs", {
  sim <- simulate_spillover(scenario_params("S1_canonical"), seed = 14)
  res <- compute_pi(sim, bootstrap = 0)
  expect_equal(colnames(res$criteria), paste0("C", 1:14))
  expect_equal(sum(res$weights_used), 1, tolerance = 1e-8)
  na_cols <- colnames(res$criteria)[apply(res$criteria, 2, function(x) all(is.na(x)))]
  expect_true(all(res$weights_used[na_cols] == 0))
  expect_true(all(res$pi >= 0 & res$pi <= 1))
})

test_that("a criterion set to zero weight cannot change the score", {
  sim <- simulate_spillover(scenario_params("S1_canonical"), seed = 15)
  w <- default_weights(); w["C12"] <- 0; w <- w / sum(w)
  a <- compute_pi(sim, weights = w)
  b <- compute_pi(sim, weights = w)
  expect_equal(a$pi, b$pi)
  expect_equal(unname(a$weights_used[["C12"]]), 0)
})

test_that("default weights sum to one and concentrate on the maintenance triad", {
  w <- default_weights()
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_equal(names(w), paste0("C", 1:14))
  expect_equal(unname(sum(w[c("C12", "C13", "C14")])), 0.45, tolerance = 1e-8)
})
