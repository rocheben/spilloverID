#' Default parameters for the simulation engine
#'
#' @return named list of default parameter values.
#' @export
scenario_params <- function(name = "S1_canonical") {
  p <- list(
    n_sites = 30, n_hosts = 4, true_reservoir = 1, bridge_host = NA,
    n_humans_per_site = 60, n_animals_per_host_site = 40,
    n_years = 15, seasons_per_year = 4,
    animal_time_resolution = 2,
    anthrop_effect_on_abundance = c(0.5, -0.3, 0.0, 0.2),
    anthrop_effect_on_exposure = 0.4,
    anthrop_effect_on_prevalence = 0.3,
    base_abundance = c(100, 80, 60, 40),
    R0_reservoir = 2.5, cross_species_beta = 0.05, recovery_rate = 0.4,
    seasonal_amplitude = 0.4, seasonal_phase = 0.0, env_noise_sd = 0.15,
    alpha_reservoir = 0.05, alpha_bridge = 0.00,
    alpha_environment = 0.003, alpha_cryptic_bg = 0.0005,
    age_min = 5, age_max = 75, residence_share_local = 0.7,
    serology_sensitivity = 0.92, serology_specificity = 0.97,
    exposure_misclass_rate = 0.10, antibody_persistence = 999
  )
  switch(name,
    "S1_canonical" = p,
    "S2_bridge" = {
      p$bridge_host <- 2; p$alpha_reservoir <- 0.035; p$alpha_bridge <- 0.025; p
    },
    "S3_cryptic" = {
      p$alpha_reservoir <- 0.015; p$alpha_environment <- 0.025
      p$alpha_cryptic_bg <- 0.003; p
    },
    "S4_confounded" = {
      p$anthrop_effect_on_abundance <- c(1.0, -0.5, 0.0, 0.3)
      p$anthrop_effect_on_exposure <- 0.8; p$anthrop_effect_on_prevalence <- 0.7; p
    },
    "S5_codistributed" = {
      p$anthrop_effect_on_abundance <- c(0.5, 0.5, 0.0, -0.2); p
    },
    "S6_strong_seasonal" = {
      p$seasonal_amplitude <- 0.7; p$env_noise_sd <- 0.25; p
    },
    "S7_low_power" = {
      p$n_humans_per_site <- 20; p$n_animals_per_host_site <- 15
      p$n_sites <- 15; p
    },
    "S8_measurement_error" = {
      p$exposure_misclass_rate <- 0.30; p$serology_sensitivity <- 0.80
      p$serology_specificity <- 0.90; p
    },
    "S9_weak_signal" = {
      p$alpha_reservoir <- 0.010; p$alpha_cryptic_bg <- 0.002; p
    },
    "S10_short_immunity" = {
      p$antibody_persistence <- 3.0; p
    },
    "S11_no_spillover" = {
      # Negative control: no host-driven transmission anywhere.
      # Human seropositivity is driven only by a low cryptic background
      # (representing assay false-positives, unrelated exposures, etc.).
      p$alpha_reservoir <- 0
      p$alpha_bridge <- 0
      p$alpha_environment <- 0
      p$alpha_cryptic_bg <- 0.0008
      p$true_reservoir <- NA  # no true reservoir
      p
    },
    stop("Unknown scenario: ", name)
  )
}

#' Simulate a multi-host wildlife-human spillover dataset
#'
#' Implements a SIRS dynamic for the reservoir with seasonal forcing and
#' stochastic noise, spillover-driven dynamics for accidental hosts and
#' optionally a bridge host, and a human population with individualized
#' exposure history.
#'
#' @param params list of parameters; see \code{\link{scenario_params}}.
#' @param seed integer for reproducibility.
#'
#' @return A \code{spillover_data} object with truth metadata in
#'   \code{attr(result, "truth")}.
#' @export
simulate_spillover <- function(params = scenario_params(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- params$n_hosts; S <- params$n_sites
  k_star <- params$true_reservoir
  k_bridge <- params$bridge_host
  no_spillover <- is.na(k_star)
  n_years <- params$n_years; sps <- params$seasons_per_year
  n_steps <- n_years * sps; dt <- 1 / sps

  A_s <- stats::runif(S, 0, 1)

  abundance <- matrix(NA, S, K)
  for (k in seq_len(K)) {
    log_a <- log(params$base_abundance[k]) +
             params$anthrop_effect_on_abundance[k] * A_s
    abundance[, k] <- pmax(stats::rpois(S, exp(log_a)), 1)
  }

  # Reservoir SIRS dynamics (skip entirely in the no-spillover case)
  I_res <- matrix(0, n_steps, S)
  if (!no_spillover) {
    p_init <- 1 - 1 / params$R0_reservoir
    I_res[1, ] <- pmin(pmax(p_init * (1 + stats::rnorm(S, 0, 0.1)), 0.01), 0.9)
    beta_base <- params$R0_reservoir * params$recovery_rate
    for (t in 2:n_steps) {
      season <- sin(2 * pi * (t / sps) + params$seasonal_phase)
      beta_t <- beta_base * (1 + params$seasonal_amplitude * season)
      beta_t_sites <- pmax(beta_t * (1 + params$anthrop_effect_on_prevalence * (A_s - 0.5)), 0.01)
      noise <- stats::rnorm(S, 0, params$env_noise_sd)
      susc <- 1 - I_res[t-1, ]
      new_inf <- beta_t_sites * susc * I_res[t-1, ] * dt * exp(noise)
      recov <- params$recovery_rate * I_res[t-1, ] * dt
      I_res[t, ] <- pmin(pmax(I_res[t-1, ] + new_inf - recov, 0.001), 0.95)
    }
  }

  I_hosts <- array(0, dim = c(n_steps, S, K))
  if (!no_spillover) I_hosts[, , k_star] <- I_res
  for (k in seq_len(K)) {
    if (!no_spillover && k == k_star) next
    for (t in 2:n_steps) {
      susc <- 1 - I_hosts[t-1, , k]
      if (!is.na(k_bridge) && k == k_bridge) {
        from_res <- 0.4 * I_res[t-1, ] * susc * dt
        self_t <- 0.5 * params$recovery_rate * I_hosts[t-1, , k] * susc * dt
        recov <- params$recovery_rate * 1.5 * I_hosts[t-1, , k] * dt
        I_hosts[t, , k] <- pmin(pmax(I_hosts[t-1, , k] + from_res + self_t - recov, 0), 0.95)
      } else {
        from_res <- params$cross_species_beta * I_res[t-1, ] * susc * dt
        recov <- params$recovery_rate * 2 * I_hosts[t-1, , k] * dt
        I_hosts[t, , k] <- pmin(pmax(I_hosts[t-1, , k] + from_res - recov, 0), 0.95)
      }
    }
  }

  # Animal seroprevalence trajectory
  sero_pop <- array(0, dim = c(n_steps, S, K))
  sero_pop[1, , ] <- 0.1 * I_hosts[1, , ]
  for (t in 2:n_steps) {
    for (k in seq_len(K)) {
      lam_k <- I_hosts[t-1, , k] * params$recovery_rate * 3
      sero_pop[t, , k] <- pmin(pmax(sero_pop[t-1, , k] +
                                    (1 - sero_pop[t-1, , k]) * lam_k * dt, 0), 0.99)
    }
  }

  # End-time sampling
  prev_pcr_true <- I_hosts[n_steps, , ]
  prev_sero_true <- sero_pop[n_steps, , ]
  # In no-spillover negative control, add a small assay-background prevalence
  # to mimic realistic cross-reactivity / specificity issues.
  if (no_spillover) {
    prev_sero_true[, ] <- pmin(prev_sero_true[, ] + 0.02, 0.10)
  }
  animal_rows <- list()
  i <- 1
  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      n <- params$n_animals_per_host_site
      n_pcr <- stats::rbinom(1, n, prev_pcr_true[s, k])
      n_sero <- stats::rbinom(1, n, prev_sero_true[s, k])
      animal_rows[[i]] <- data.frame(site = s, host = k, n_sampled = n,
                                     n_pcr = n_pcr, prev_pcr = n_pcr / n,
                                     n_sero = n_sero, prev_sero = n_sero / n)
      i <- i + 1
    }
  }
  animal_df <- do.call(rbind, animal_rows)

  # Animal individual data
  animal_indiv_rows <- list()
  i <- 1
  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      n <- params$n_animals_per_host_site
      age <- stats::runif(n, 0.1, 2.5)
      if (!no_spillover && k == k_star) {
        foi_k <- -log(1 - prev_sero_true[s, k]) / 1
      } else {
        foi_k <- -log(1 - prev_sero_true[s, k] + 1e-4) / 1 * 0.5
      }
      p_indiv <- 1 - exp(-foi_k * age)
      sero_true <- stats::rbinom(n, 1, p_indiv)
      sero_obs <- ifelse(sero_true == 1,
                         stats::rbinom(n, 1, params$serology_sensitivity),
                         stats::rbinom(n, 1, 1 - params$serology_specificity))
      animal_indiv_rows[[i]] <- data.frame(site = s, host = k, age = age,
                                           sero = sero_obs, sero_true = sero_true)
      i <- i + 1
    }
  }
  animal_indiv_df <- do.call(rbind, animal_indiv_rows)

  # Animal time series
  obs_per_year <- params$animal_time_resolution
  obs_step <- max(1, sps %/% obs_per_year)
  obs_times <- seq(1, n_steps, by = obs_step)
  animal_ts_rows <- list(); i <- 1
  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      for (t in obs_times) {
        animal_ts_rows[[i]] <- data.frame(site = s, host = k,
                                          year = (t - 1) / sps,
                                          prev_pcr_true = I_hosts[t, s, k])
        i <- i + 1
      }
    }
  }
  animal_ts <- do.call(rbind, animal_ts_rows)

  # Human population
  human_rows <- list(); ind_idx <- 1
  for (s in seq_len(S)) {
    n <- params$n_humans_per_site
    age <- stats::runif(n, params$age_min, params$age_max)
    sex <- stats::rbinom(n, 1, 0.5)
    prof <- stats::rbinom(n, 1, 0.3)
    residence <- pmin(age, pmax(age * params$residence_share_local +
                                stats::rnorm(n, 0, 3), 1))
    exp_env_cont <- stats::plogis(-1 + 1.5 * A_s[s] + stats::rnorm(n, 0, 0.3))
    exposures_true <- matrix(0, n, K)
    for (k in seq_len(K)) {
      logit_e <- -2 + 0.3 * log(abundance[s, k] + 1) + 1.2 * prof + 0.5 * A_s[s] +
                 stats::rnorm(n, 0, 0.5)
      exposures_true[, k] <- stats::plogis(logit_e)
    }
    household <- sample.int(max(1, n %/% 4), n, replace = TRUE)
    for (i in seq_len(n)) {
      entry_year <- n_years - residence[i]
      entry_step <- max(1, floor(entry_year * sps))
      sero_true <- 0; inf_time <- NA_real_
      for (t in entry_step:n_steps) {
        if (no_spillover) {
          lam_t <- params$alpha_cryptic_bg
        } else {
          lam_t <- params$alpha_reservoir * exposures_true[i, k_star] * I_hosts[t, s, k_star] +
                   params$alpha_environment * exp_env_cont[i] * I_hosts[t, s, k_star] +
                   params$alpha_cryptic_bg
        }
        if (!is.na(k_bridge)) {
          lam_t <- lam_t + params$alpha_bridge * exposures_true[i, k_bridge] *
                          I_hosts[t, s, k_bridge]
        }
        p_step <- 1 - exp(-lam_t * dt)
        if (stats::runif(1) < p_step) {
          sero_true <- 1; inf_time <- (t - 1) / sps
          if (n_years - inf_time > params$antibody_persistence) sero_true <- 0
          break
        }
      }
      sero_obs <- if (sero_true == 1) stats::rbinom(1, 1, params$serology_sensitivity)
                  else stats::rbinom(1, 1, 1 - params$serology_specificity)
      row <- data.frame(site = s, ind = ind_idx, age = age[i], sex = sex[i],
                        prof = prof[i], residence = residence[i],
                        household = household[i], anthrop = A_s[s],
                        env_cont = exp_env_cont[i],
                        env_declared = as.integer(exp_env_cont[i] > 0.5),
                        sero = sero_obs, sero_true = sero_true,
                        inf_time = inf_time)
      for (k in seq_len(K)) {
        base <- as.integer(exposures_true[i, k] > 0.3)
        if (stats::runif(1) < params$exposure_misclass_rate) base <- 1 - base
        row[[paste0("exp_h", k)]] <- base
      }
      human_rows[[ind_idx]] <- row
      ind_idx <- ind_idx + 1
    }
  }
  human_df <- do.call(rbind, human_rows)

  # Human incidence time series
  inc_rows <- list(); i <- 1
  for (s in seq_len(S)) {
    sub <- human_df[human_df$site == s, ]
    for (y in 0:(n_years - 1)) {
      n_inf <- sum(!is.na(sub$inf_time) & sub$inf_time >= y & sub$inf_time < y + 1)
      inc_rows[[i]] <- data.frame(site = s, year = y, n_new_inf = n_inf,
                                  incidence = n_inf / max(nrow(sub), 1))
      i <- i + 1
    }
  }
  human_ts <- do.call(rbind, inc_rows)

  # Site df
  site_df <- data.frame(site = seq_len(S), anthrop = A_s)
  for (k in seq_len(K)) {
    site_df[[paste0("prev_pcr_h", k)]] <- animal_df$prev_pcr[animal_df$host == k]
    site_df[[paste0("prev_sero_h", k)]] <- animal_df$prev_sero[animal_df$host == k]
    site_df[[paste0("abund_h", k)]] <- abundance[, k]
  }

  obj <- spillover_data(human = human_df, animal = animal_df,
                        animal_indiv = animal_indiv_df,
                        animal_ts = animal_ts, human_ts = human_ts,
                        site = site_df, n_hosts = K,
                        host_names = paste0("h", seq_len(K)))

  # Derive animal_periods (epidemic vs inter-epidemic) from time series.
  # An "epidemic year" is one where human incidence exceeds the median across years.
  inc_per_year <- aggregate(incidence ~ year, data = human_ts, FUN = mean)
  median_inc <- stats::median(inc_per_year$incidence, na.rm = TRUE)
  epi_years <- inc_per_year$year[inc_per_year$incidence > median_inc]
  inter_years <- inc_per_year$year[inc_per_year$incidence <= median_inc]
  ap_rows <- list()
  for (k in seq_len(K)) {
    epi_rows  <- animal_ts[animal_ts$host == k & animal_ts$year %in% epi_years, ]
    inter_rows <- animal_ts[animal_ts$host == k & animal_ts$year %in% inter_years, ]
    if (nrow(epi_rows) > 0) {
      p_epi <- mean(epi_rows$prev_pcr_true, na.rm = TRUE)
      n_per_epi <- 60
      ap_rows[[length(ap_rows) + 1]] <- data.frame(
        site = 1, host = k, period = "epidemic",
        n_sampled = n_per_epi,
        n_sero = round(min(max(p_epi, 0), 1) * n_per_epi)
      )
    }
    if (nrow(inter_rows) > 0) {
      p_inter <- mean(inter_rows$prev_pcr_true, na.rm = TRUE)
      n_per_inter <- 60
      ap_rows[[length(ap_rows) + 1]] <- data.frame(
        site = 1, host = k, period = "interepidemic",
        n_sampled = n_per_inter,
        n_sero = round(min(max(p_inter, 0), 1) * n_per_inter)
      )
    }
  }
  if (length(ap_rows) > 0) obj$animal_periods <- do.call(rbind, ap_rows)

  attr(obj, "truth") <- list(k_star = k_star, k_bridge = k_bridge,
                             prevalence_pcr = prev_pcr_true,
                             prevalence_sero = prev_sero_true,
                             A_s = A_s, abundance = abundance)
  obj
}
