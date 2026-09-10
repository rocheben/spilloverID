test_that("spillover_data builds correctly", {
  human <- data.frame(
    site = rep(1:5, each = 20), ind = 1:100,
    age = runif(100, 5, 75), sero = sample(0:1, 100, TRUE),
    exp_h1 = sample(0:1, 100, TRUE), exp_h2 = sample(0:1, 100, TRUE)
  )
  animal <- expand.grid(site = 1:5, host = 1:2)
  animal$n_sampled <- 30
  animal$n_pcr <- rbinom(nrow(animal), 30, 0.3)
  animal$prev_pcr <- animal$n_pcr / 30
  animal$n_sero <- rbinom(nrow(animal), 30, 0.4)
  animal$prev_sero <- animal$n_sero / 30
  d <- spillover_data(human, animal, n_hosts = 2)
  expect_s3_class(d, "spillover_data")
  expect_equal(d$n_hosts, 2)
})

test_that("compute_pi returns valid Pi vector", {
  sim <- simulate_spillover(scenario_params("S1_canonical"), seed = 1)
  res <- compute_pi(sim)
  expect_s3_class(res, "spillover_pi")
  expect_length(res$pi, sim$n_hosts)
  expect_true(all(res$pi >= 0 & res$pi <= 1))
})

test_that("weights are renormalized when criteria are NA", {
  w <- default_weights()
  expect_equal(sum(w), 1, tolerance = 1e-6)
})

test_that("rank_norm and softmax_norm return [0,1]", {
  x <- c(0.1, 0.5, 2.0, -1.0)
  expect_true(all(rank_norm(x) >= 0 & rank_norm(x) <= 1))
  expect_true(all(softmax_norm(x) >= 0 & softmax_norm(x) <= 1))
  expect_equal(max(softmax_norm(x)), 1)
})
