# =====================================================================
# simulate_ext.R -- Additive extension of the spilloverID simulator.
#
# Adds the scenario families that relax the single-species assumption:
#   S12_congeners      two co-located congeneric hosts, near-identical
#                      maintenance dynamics, serologically cross-reactive
#   S13_dual_symmetric two hosts genuinely co-maintain the pathogen (50/50)
#   S13b_dual_asym     two co-maintaining hosts with asymmetric spillover
#   S14_hidden_host    the true reservoir is NOT among the sampled hosts
#   S1c_ref / S11c_ref reference arms re-run through THIS generator so that
#                      the new scenarios are compared to a like-for-like
#                      baseline rather than to the published run
#
# DESIGN NOTE. The published simulate_spillover() hard-codes a single
# reservoir index (params$true_reservoir, scalar) and has no notion of a
# taxonomic grouping, of assay cross-reactivity between candidate hosts, or
# of unsampled hosts. Those are structural, not parametric, limitations, so
# this file provides a *generalised* generator, simulate_multihost(), rather
# than a wrapper. It is a deliberate fork of simulate_spillover(): the
# dynamics, sampling, observation model and animal_periods construction are
# reproduced line-for-line except where a generalisation is required, and
# validate_ext.R checks that the fork reproduces the published S1/S11
# behaviour when run with legacy settings. NO existing file is modified.
#
# Generalisations, each switched off by default:
#   maintenance   data.frame(host, R0, share): every host that genuinely
#                 maintains the pathogen. `share` is that host's share of
#                 the human force of infection (share = 0 => the host
#                 maintains the agent but does not seed human cases).
#   genus         integer vector of length n_hosts, taxonomic grouping.
#   crossreact    list(pair = c(a, b), xi = numeric): the serological assay
#                 misassigns a fraction xi of one congener's signal to the
#                 other, symmetrically. PCR is left species-specific.
#   exposure_confusion  probability that a respondent cannot attribute
#                 contact to one member of the congener pair (both
#                 exposure indicators are then set to their union).
#   hidden_hosts  hosts simulated but deleted from the observed data set
#                 before scoring (remaining hosts are renumbered).
#   abundance_shared_latent  site-level latent shared by a pair of hosts,
#                 i.e. genuine co-location.
# =====================================================================

# ---------------------------------------------------------------------
# Scenario parameters. Legacy names are delegated untouched to the
# published scenario_params(); only new names are handled here.
# ---------------------------------------------------------------------
scenario_params_ext <- function(name = "S12_congeners") {

  legacy <- c("S1_canonical", "S2_bridge", "S3_cryptic", "S4_confounded",
              "S5_codistributed", "S6_strong_seasonal", "S7_low_power",
              "S8_measurement_error", "S9_weak_signal", "S10_short_immunity",
              "S11_no_spillover")
  if (name %in% legacy) {
    p <- scenario_params(name)
    # Express the legacy scalar truth in the generalised vocabulary.
    p$maintenance <- if (is.na(p$true_reservoir)) {
      data.frame(host = integer(0), R0 = double(0), share = double(0))
    } else {
      data.frame(host = p$true_reservoir, R0 = p$R0_reservoir, share = 1)
    }
    p$genus <- seq_len(p$n_hosts)
    p$crossreact <- NULL
    p$exposure_confusion <- 0
    p$hidden_hosts <- integer(0)
    p$abundance_shared_latent <- NULL
    p$congener_similarity_range <- NULL
    p$share_range <- NULL
    return(p)
  }

  # Base template: identical to the published defaults, extended to K = 5.
  p <- list(
    n_sites = 30, n_hosts = 5, true_reservoir = 1, bridge_host = NA,
    n_humans_per_site = 60, n_animals_per_host_site = 40,
    n_years = 15, seasons_per_year = 4,
    animal_time_resolution = 2,
    anthrop_effect_on_abundance = c(0.5, 0.5, -0.3, 0.0, 0.2),
    anthrop_effect_on_exposure = 0.4,
    anthrop_effect_on_prevalence = 0.3,
    base_abundance = c(100, 95, 80, 60, 40),
    R0_reservoir = 2.5, cross_species_beta = 0.05, recovery_rate = 0.4,
    seasonal_amplitude = 0.4, seasonal_phase = 0.0, env_noise_sd = 0.15,
    alpha_reservoir = 0.05, alpha_bridge = 0.00,
    alpha_environment = 0.003, alpha_cryptic_bg = 0.0005,
    age_min = 5, age_max = 75, residence_share_local = 0.7,
    serology_sensitivity = 0.92, serology_specificity = 0.97,
    exposure_misclass_rate = 0.10, antibody_persistence = 999,
    # --- generalised fields, inert unless set ---
    maintenance = data.frame(host = 1, R0 = 2.5, share = 1),
    genus = 1:5,
    crossreact = NULL,
    exposure_confusion = 0,
    hidden_hosts = integer(0),
    abundance_shared_latent = NULL,
    congener_similarity_range = NULL,
    share_range = NULL,
    crossreact_applies_to_periods = TRUE
  )

  switch(name,

    # ---- S12: two congeners of the same genus -----------------------
    # h1 and h2 are congeneric (genus 1), co-located, and both maintain
    # the agent; h2's maintenance intensity is rho x that of h1 with
    # rho drawn per replicate in [0.70, 0.90]. Only h1 seeds human cases
    # (share = c(1, 0)) -- that is the ground truth of interest: which
    # *species* of the genus is responsible. The serological assay does
    # not resolve the pair (xi), nor do respondents (exposure_confusion).
    "S12_congeners" = {
      p$genus <- c(1, 1, 2, 3, 4)
      p$maintenance <- data.frame(host = c(1, 2), R0 = c(2.5, NA),
                                  share = c(1, 0))
      p$congener_similarity_range <- c(0.70, 0.90)   # rho, drawn per replicate
      p$crossreact <- list(pair = c(1, 2), xi = 0.35)
      p$exposure_confusion <- 0.30
      p$abundance_shared_latent <- list(pair = c(1, 2), sd = 0.30)
      p$base_abundance <- c(100, 95, 80, 60, 40)
      p$anthrop_effect_on_abundance <- c(0.5, 0.5, -0.3, 0.0, 0.2)
      p
    },

    # ---- S12b: same congener ecology, but a species-resolving assay --
    # Identical to S12 in every respect except that the serological assay
    # discriminates the two congeners (xi = 0) and respondents attribute
    # contact correctly (tau = 0). Contrasting S12b with S12 partitions the
    # loss of species resolution into (i) irreducible ecological similarity
    # and (ii) the measurement layer -- which is the actionable half.
    "S12b_congeners_resolved" = {
      p$genus <- c(1, 1, 2, 3, 4)
      p$maintenance <- data.frame(host = c(1, 2), R0 = c(2.5, NA),
                                  share = c(1, 0))
      p$congener_similarity_range <- c(0.70, 0.90)
      p$crossreact <- list(pair = c(1, 2), xi = 0.0)
      p$exposure_confusion <- 0
      p$abundance_shared_latent <- list(pair = c(1, 2), sd = 0.30)
      p$base_abundance <- c(100, 95, 80, 60, 40)
      p$anthrop_effect_on_abundance <- c(0.5, 0.5, -0.3, 0.0, 0.2)
      p
    },

    # ---- S13: two genuine, symmetric co-reservoirs -------------------
    # h1 and h2 are NOT congeners (different genera) and both truly
    # maintain the pathogen and both seed human cases, equally.
    # There is no unique k_star.
    "S13_dual_symmetric" = {
      p$genus <- c(1, 2, 3, 4, 5)
      p$maintenance <- data.frame(host = c(1, 2), R0 = c(2.5, 2.5),
                                  share = c(0.5, 0.5))
      p$base_abundance <- c(100, 90, 80, 60, 40)
      p$anthrop_effect_on_abundance <- c(0.5, 0.4, -0.3, 0.0, 0.2)
      p
    },

    # ---- S13b: two co-reservoirs, asymmetric spillover ---------------
    # Same as S13 but h1 contributes a fraction s ~ U(0.60, 0.90) of the
    # human force of infection. Gives a continuum from "co-equal" to
    # "one dominant reservoir + one secondary".
    "S13b_dual_asymmetric" = {
      p$genus <- c(1, 2, 3, 4, 5)
      p$maintenance <- data.frame(host = c(1, 2), R0 = c(2.5, 2.5),
                                  share = c(NA, NA))
      p$share_range <- c(0.60, 0.90)     # share of h1, drawn per replicate
      p$base_abundance <- c(100, 90, 80, 60, 40)
      p$anthrop_effect_on_abundance <- c(0.5, 0.4, -0.3, 0.0, 0.2)
      p
    },

    # ---- S14: incomplete inventory (Drexler's squirrel) --------------
    # Five hosts are simulated; host 5 is the sole true reservoir and is
    # deleted from the data before scoring. The four hosts that ARE
    # sampled are seropositive only through cross-species spillover from
    # the missing reservoir. A well-behaved score should look like the
    # negative control here.
    "S14_hidden_reservoir" = {
      p$genus <- c(1, 2, 3, 4, 5)
      p$maintenance <- data.frame(host = 5, R0 = 2.5, share = 1)
      p$hidden_hosts <- 5L
      p$base_abundance <- c(80, 60, 50, 40, 100)
      p$anthrop_effect_on_abundance <- c(-0.3, 0.0, 0.2, 0.1, 0.5)
      p
    },

    # ---- reference arms, run through THIS generator ------------------
    "S1c_ref_single" = {
      p$n_hosts <- 4
      p$genus <- 1:4
      p$maintenance <- data.frame(host = 1, R0 = 2.5, share = 1)
      p$base_abundance <- c(100, 80, 60, 40)
      p$anthrop_effect_on_abundance <- c(0.5, -0.3, 0.0, 0.2)
      p
    },
    "S11c_ref_null" = {
      p$n_hosts <- 4
      p$genus <- 1:4
      p$maintenance <- data.frame(host = integer(0), R0 = double(0),
                                  share = double(0))
      p$base_abundance <- c(100, 80, 60, 40)
      p$anthrop_effect_on_abundance <- c(0.5, -0.3, 0.0, 0.2)
      p$alpha_reservoir <- 0
      p$alpha_bridge <- 0
      p$alpha_environment <- 0
      p$alpha_cryptic_bg <- 0.0008
      p$true_reservoir <- NA
      p
    },
    stop("Unknown extended scenario: ", name)
  )
}


# ---------------------------------------------------------------------
# Generalised generator.
# ---------------------------------------------------------------------
simulate_multihost <- function(params, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  K <- params$n_hosts; S <- params$n_sites
  n_years <- params$n_years; sps <- params$seasons_per_year
  n_steps <- n_years * sps; dt <- 1 / sps
  k_bridge <- params$bridge_host

  maint <- params$maintenance
  # Per-replicate draws (must happen after set.seed so they are reproducible)
  rho <- NA_real_
  if (!is.null(params$congener_similarity_range)) {
    rho <- stats::runif(1, params$congener_similarity_range[1],
                           params$congener_similarity_range[2])
    # congener's maintenance intensity is rho x the focal host's
    maint$R0[2] <- 1 + rho * (maint$R0[1] - 1)
  }
  if (!is.null(params$share_range)) {
    s1 <- stats::runif(1, params$share_range[1], params$share_range[2])
    maint$share <- c(s1, 1 - s1)
  }
  no_spillover <- nrow(maint) == 0 || sum(maint$share) == 0
  maint_hosts <- maint$host

  A_s <- stats::runif(S, 0, 1)

  # ---- host abundance, with an optional shared site latent for a pair
  shared_eps <- rep(0, S)
  pair_lat <- params$abundance_shared_latent
  if (!is.null(pair_lat)) shared_eps <- stats::rnorm(S, 0, pair_lat$sd)
  abundance <- matrix(NA, S, K)
  for (k in seq_len(K)) {
    log_a <- log(params$base_abundance[k]) +
             params$anthrop_effect_on_abundance[k] * A_s
    if (!is.null(pair_lat) && k %in% pair_lat$pair) log_a <- log_a + shared_eps
    abundance[, k] <- pmax(stats::rpois(S, exp(log_a)), 1)
  }

  # ---- SIRS dynamics, one trajectory per maintenance host
  I_hosts <- array(0, dim = c(n_steps, S, K))
  for (j in seq_len(nrow(maint))) {
    k <- maint$host[j]; R0k <- maint$R0[j]
    p_init <- 1 - 1 / R0k
    I_hosts[1, , k] <- pmin(pmax(p_init * (1 + stats::rnorm(S, 0, 0.1)), 0.01), 0.9)
    beta_base <- R0k * params$recovery_rate
    for (t in 2:n_steps) {
      season <- sin(2 * pi * (t / sps) + params$seasonal_phase)
      beta_t <- beta_base * (1 + params$seasonal_amplitude * season)
      beta_t_sites <- pmax(beta_t * (1 + params$anthrop_effect_on_prevalence *
                                       (A_s - 0.5)), 0.01)
      noise <- stats::rnorm(S, 0, params$env_noise_sd)
      susc <- 1 - I_hosts[t-1, , k]
      new_inf <- beta_t_sites * susc * I_hosts[t-1, , k] * dt * exp(noise)
      recov <- params$recovery_rate * I_hosts[t-1, , k] * dt
      I_hosts[t, , k] <- pmin(pmax(I_hosts[t-1, , k] + new_inf - recov, 0.001), 0.95)
    }
  }

  # Aggregate reservoir pressure felt by non-maintenance hosts: the mean
  # over maintenance hosts (so that adding a second reservoir does not
  # mechanically double the prevalence of accidental hosts, which would
  # confound the multi-reservoir question with a pure dose effect).
  I_agg <- matrix(0, n_steps, S)
  if (length(maint_hosts) > 0) {
    for (k in maint_hosts) I_agg <- I_agg + I_hosts[, , k]
    I_agg <- I_agg / length(maint_hosts)
  }

  for (k in seq_len(K)) {
    if (k %in% maint_hosts) next
    for (t in 2:n_steps) {
      susc <- 1 - I_hosts[t-1, , k]
      if (!is.na(k_bridge) && k == k_bridge) {
        from_res <- 0.4 * I_agg[t-1, ] * susc * dt
        self_t <- 0.5 * params$recovery_rate * I_hosts[t-1, , k] * susc * dt
        recov <- params$recovery_rate * 1.5 * I_hosts[t-1, , k] * dt
        I_hosts[t, , k] <- pmin(pmax(I_hosts[t-1, , k] + from_res + self_t - recov, 0), 0.95)
      } else {
        from_res <- params$cross_species_beta * I_agg[t-1, ] * susc * dt
        recov <- params$recovery_rate * 2 * I_hosts[t-1, , k] * dt
        I_hosts[t, , k] <- pmin(pmax(I_hosts[t-1, , k] + from_res - recov, 0), 0.95)
      }
    }
  }

  # ---- animal seroprevalence trajectory
  sero_pop <- array(0, dim = c(n_steps, S, K))
  sero_pop[1, , ] <- 0.1 * I_hosts[1, , ]
  for (t in 2:n_steps) {
    for (k in seq_len(K)) {
      lam_k <- I_hosts[t-1, , k] * params$recovery_rate * 3
      sero_pop[t, , k] <- pmin(pmax(sero_pop[t-1, , k] +
                                    (1 - sero_pop[t-1, , k]) * lam_k * dt, 0), 0.99)
    }
  }

  prev_pcr_true  <- I_hosts[n_steps, , , drop = TRUE]
  prev_sero_true <- sero_pop[n_steps, , , drop = TRUE]
  if (K == 1) {
    prev_pcr_true  <- matrix(prev_pcr_true, ncol = 1)
    prev_sero_true <- matrix(prev_sero_true, ncol = 1)
  }
  if (no_spillover) prev_sero_true[, ] <- pmin(prev_sero_true[, ] + 0.02, 0.10)

  # ---- serological cross-reactivity between congeners.
  # The assay reports, for each member of the pair, a mixture of the two
  # species' true seroprevalences. PCR is left species-specific.
  xr <- params$crossreact
  prev_sero_obs <- prev_sero_true
  if (!is.null(xr)) {
    a <- xr$pair[1]; b <- xr$pair[2]; xi <- xr$xi
    prev_sero_obs[, a] <- (1 - xi) * prev_sero_true[, a] + xi * prev_sero_true[, b]
    prev_sero_obs[, b] <- (1 - xi) * prev_sero_true[, b] + xi * prev_sero_true[, a]
  }

  # ---- site x host aggregate animal sampling
  animal_rows <- list(); i <- 1
  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      n <- params$n_animals_per_host_site
      n_pcr  <- stats::rbinom(1, n, prev_pcr_true[s, k])
      n_sero <- stats::rbinom(1, n, prev_sero_obs[s, k])
      animal_rows[[i]] <- data.frame(site = s, host = k, n_sampled = n,
                                     n_pcr = n_pcr, prev_pcr = n_pcr / n,
                                     n_sero = n_sero, prev_sero = n_sero / n)
      i <- i + 1
    }
  }
  animal_df <- do.call(rbind, animal_rows)

  # ---- individual animal data (age-structured), with cross-reactivity
  animal_indiv_rows <- list(); i <- 1
  for (s in seq_len(S)) {
    for (k in seq_len(K)) {
      n <- params$n_animals_per_host_site
      age <- stats::runif(n, 0.1, 2.5)
      foi_of <- function(kk) {
        if (kk %in% maint_hosts) -log(1 - prev_sero_true[s, kk]) / 1
        else -log(1 - prev_sero_true[s, kk] + 1e-4) / 1 * 0.5
      }
      if (!is.null(xr) && k %in% xr$pair) {
        other <- setdiff(xr$pair, k)
        foi_k <- (1 - xr$xi) * foi_of(k) + xr$xi * foi_of(other)
      } else {
        foi_k <- foi_of(k)
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

  # ---- animal time series
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

  # ---- human population
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
          lam_res <- 0
          for (j in seq_len(nrow(maint))) {
            kj <- maint$host[j]
            lam_res <- lam_res + params$alpha_reservoir * maint$share[j] *
                       exposures_true[i, kj] * I_hosts[t, s, kj]
          }
          lam_t <- lam_res +
                   params$alpha_environment * exp_env_cont[i] * I_agg[t, s] +
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
      decl <- integer(K)
      for (k in seq_len(K)) {
        base <- as.integer(exposures_true[i, k] > 0.3)
        if (stats::runif(1) < params$exposure_misclass_rate) base <- 1 - base
        decl[k] <- base
      }
      # Congener exposure confusion: with probability tau the respondent
      # cannot attribute the contact to one species of the pair, and both
      # indicators take the union value.
      if (!is.null(xr) && params$exposure_confusion > 0) {
        if (stats::runif(1) < params$exposure_confusion) {
          u <- max(decl[xr$pair[1]], decl[xr$pair[2]])
          decl[xr$pair[1]] <- u; decl[xr$pair[2]] <- u
        }
      }
      for (k in seq_len(K)) row[[paste0("exp_h", k)]] <- decl[k]
      human_rows[[ind_idx]] <- row
      ind_idx <- ind_idx + 1
    }
  }
  human_df <- do.call(rbind, human_rows)

  # ---- human incidence time series
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

  # ---- animal_periods (support of C12), built exactly as in the
  # published simulator, then subjected to the same serological
  # cross-reactivity as every other serological measurement.
  inc_per_year <- aggregate(incidence ~ year, data = human_ts, FUN = mean)
  median_inc <- stats::median(inc_per_year$incidence, na.rm = TRUE)
  epi_years <- inc_per_year$year[inc_per_year$incidence > median_inc]
  inter_years <- inc_per_year$year[inc_per_year$incidence <= median_inc]
  p_epi_k <- rep(NA_real_, K); p_int_k <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    er <- animal_ts[animal_ts$host == k & animal_ts$year %in% epi_years, ]
    ir <- animal_ts[animal_ts$host == k & animal_ts$year %in% inter_years, ]
    if (nrow(er) > 0) p_epi_k[k] <- mean(er$prev_pcr_true, na.rm = TRUE)
    if (nrow(ir) > 0) p_int_k[k] <- mean(ir$prev_pcr_true, na.rm = TRUE)
  }
  if (!is.null(xr) && isTRUE(params$crossreact_applies_to_periods)) {
    a <- xr$pair[1]; b <- xr$pair[2]; xi <- xr$xi
    pe <- p_epi_k; pi_ <- p_int_k
    p_epi_k[a] <- (1 - xi) * pe[a] + xi * pe[b]
    p_epi_k[b] <- (1 - xi) * pe[b] + xi * pe[a]
    p_int_k[a] <- (1 - xi) * pi_[a] + xi * pi_[b]
    p_int_k[b] <- (1 - xi) * pi_[b] + xi * pi_[a]
  }
  ap_rows <- list()
  for (k in seq_len(K)) {
    if (!is.na(p_epi_k[k])) ap_rows[[length(ap_rows) + 1]] <- data.frame(
      site = 1, host = k, period = "epidemic", n_sampled = 60,
      n_sero = round(min(max(p_epi_k[k], 0), 1) * 60))
    if (!is.na(p_int_k[k])) ap_rows[[length(ap_rows) + 1]] <- data.frame(
      site = 1, host = k, period = "interepidemic", n_sampled = 60,
      n_sero = round(min(max(p_int_k[k], 0), 1) * 60))
  }
  animal_periods <- if (length(ap_rows) > 0) do.call(rbind, ap_rows) else NULL

  # ---- delete unsampled hosts and renumber (incomplete inventory)
  hidden <- params$hidden_hosts
  keep <- setdiff(seq_len(K), hidden)
  remap <- rep(NA_integer_, K); remap[keep] <- seq_along(keep)
  K_obs <- length(keep)
  if (length(hidden) > 0) {
    animal_df <- animal_df[animal_df$host %in% keep, ]
    animal_df$host <- remap[animal_df$host]
    animal_indiv_df <- animal_indiv_df[animal_indiv_df$host %in% keep, ]
    animal_indiv_df$host <- remap[animal_indiv_df$host]
    animal_ts <- animal_ts[animal_ts$host %in% keep, ]
    animal_ts$host <- remap[animal_ts$host]
    if (!is.null(animal_periods)) {
      animal_periods <- animal_periods[animal_periods$host %in% keep, ]
      animal_periods$host <- remap[animal_periods$host]
    }
    drop_cols <- paste0("exp_h", hidden)
    human_df <- human_df[, setdiff(names(human_df), drop_cols), drop = FALSE]
    for (k in keep) {
      if (remap[k] != k) names(human_df)[names(human_df) == paste0("exp_h", k)] <-
        paste0("exp_h_new", remap[k])
    }
    nm <- names(human_df)
    names(human_df) <- sub("^exp_h_new", "exp_h", nm)
  }

  # ---- site-level frame (built from the OBSERVED hosts only)
  site_df <- data.frame(site = seq_len(S), anthrop = A_s)
  for (kk in seq_len(K_obs)) {
    site_df[[paste0("prev_pcr_h", kk)]]  <- animal_df$prev_pcr[animal_df$host == kk]
    site_df[[paste0("prev_sero_h", kk)]] <- animal_df$prev_sero[animal_df$host == kk]
    site_df[[paste0("abund_h", kk)]]     <- abundance[, keep[kk]]
  }

  obj <- spillover_data(human = human_df, animal = animal_df,
                        animal_indiv = animal_indiv_df,
                        animal_ts = animal_ts, human_ts = human_ts,
                        site = site_df, n_hosts = K_obs,
                        host_names = paste0("h", seq_len(K_obs)))
  if (!is.null(animal_periods)) obj$animal_periods <- animal_periods

  # ---- ground truth, expressed for a possibly non-scalar reservoir set
  k_star_obs <- remap[maint$host[maint$share > 0]]
  k_star_obs <- k_star_obs[!is.na(k_star_obs)]
  maint_obs  <- remap[maint$host]; maint_obs <- maint_obs[!is.na(maint_obs)]
  attr(obj, "truth") <- list(
    k_star_set        = k_star_obs,            # hosts that seed human cases AND are observed
    k_star_primary    = if (length(k_star_obs)) k_star_obs[which.max(maint$share[maint$share > 0])] else NA_integer_,
    maintenance_obs   = maint_obs,             # hosts that maintain AND are observed
    shares            = maint$share,
    R0                = maint$R0,
    genus_obs         = params$genus[keep],
    hidden_hosts      = hidden,
    reservoir_hidden  = length(setdiff(maint$host[maint$share > 0], keep)) > 0,
    rho_congener      = rho,
    n_hosts_observed  = K_obs,
    prevalence_pcr    = prev_pcr_true,
    prevalence_sero   = prev_sero_true,
    A_s = A_s, abundance = abundance
  )
  obj
}
