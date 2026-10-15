#' ============================================================================
#' Reservoir identification on three emblematic zoonotic outbreaks
#' with the 14-criterion composite plausibility score (spilloverID).
#' Reproducibility: set.seed(20260524); R 4.3+; spilloverID >= 0.2.0.
#' ============================================================================

suppressPackageStartupMessages({
  library(spilloverID)
  library(dplyr)
  library(jsonlite)
})

set.seed(20260524)

# ---- Dataset constructor -----------------------------------------------------
build_outbreak_dataset <- function(
  n_sites, host_names, host_prev_pcr, host_prev_sero, host_n_sampled,
  human_per_site, human_exposure_probs, spillover_intensity, exposure_rr,
  interepi_prev_ratio, endemic_signature_str, anthrop_vec = NULL
) {
  K <- length(host_names)
  if (is.null(anthrop_vec)) anthrop_vec <- runif(n_sites, 0.1, 0.9)
  animal_rows <- list()
  for (k in seq_len(K)) for (s in seq_len(n_sites)) {
    animal_rows[[length(animal_rows) + 1]] <- data.frame(
      site = s, host = k,
      n_sampled = host_n_sampled[[k]][s],
      n_pcr = round(host_prev_pcr[[k]][s] * host_n_sampled[[k]][s]),
      prev_pcr = host_prev_pcr[[k]][s],
      n_sero = round(host_prev_sero[[k]][s] * host_n_sampled[[k]][s]),
      prev_sero = host_prev_sero[[k]][s]
    )
  }
  animal <- do.call(rbind, animal_rows)

  animal_indiv_rows <- list()
  for (k in seq_len(K)) for (s in seq_len(n_sites)) {
    n <- max(30, round(host_n_sampled[[k]][s]))
    ages <- runif(n, 0.2, 5)
    str_end <- endemic_signature_str[[k]]
    base_p <- host_prev_sero[[k]][s]
    lambda_endemic <- -log(1 - min(base_p, 0.9)) / 2.5
    p_endemic <- 1 - exp(-lambda_endemic * ages)
    p_epidemic <- rep(base_p, n)
    p_age <- str_end * p_endemic + (1 - str_end) * p_epidemic
    sero <- rbinom(n, 1, pmin(p_age, 1))
    animal_indiv_rows[[length(animal_indiv_rows) + 1]] <- data.frame(
      site = s, host = k, age = ages, sero = sero
    )
  }
  animal_indiv <- do.call(rbind, animal_indiv_rows)

  human_rows <- list()
  ind_counter <- 0
  for (s in seq_len(n_sites)) for (i in seq_len(human_per_site)) {
    ind_counter <- ind_counter + 1
    age <- sample(18:75, 1)
    sex <- sample(c("M", "F"), 1, prob = c(0.55, 0.45))
    residence <- min(age - 17, sample(1:40, 1))
    exps <- sapply(seq_len(K), function(k) rbinom(1, 1, human_exposure_probs[[k]][s]))
    hazard <- spillover_intensity
    for (k in seq_len(K)) if (exps[k] == 1) hazard <- hazard * exposure_rr[[k]][s]
    sero <- rbinom(1, 1, 1 - exp(-hazard * residence))
    row <- data.frame(
      site = s, ind = ind_counter, age = age, sex = sex,
      residence = residence, sero = sero, anthrop = anthrop_vec[s],
      prof = sample(0:1, 1, prob = c(0.6, 0.4)),
      env_declared = rbinom(1, 1, 0.3)
    )
    for (k in seq_len(K)) row[[paste0("exp_h", k)]] <- exps[k]
    human_rows[[length(human_rows) + 1]] <- row
  }
  human <- do.call(rbind, human_rows)

  site_df <- data.frame(site = seq_len(n_sites), anthrop = anthrop_vec)
  for (k in seq_len(K)) {
    site_df[[paste0("prev_pcr_h", k)]] <- host_prev_pcr[[k]]
    site_df[[paste0("prev_sero_h", k)]] <- host_prev_sero[[k]]
  }

  ap_rows <- list()
  for (k in seq_len(K)) {
    p_epi <- mean(host_prev_sero[[k]])
    p_inter <- p_epi * interepi_prev_ratio[[k]]
    n_per <- 80
    ap_rows[[length(ap_rows) + 1]] <- data.frame(
      site = 1, host = k, period = "epidemic",
      n_sampled = n_per, n_sero = round(p_epi * n_per)
    )
    ap_rows[[length(ap_rows) + 1]] <- data.frame(
      site = 1, host = k, period = "interepidemic",
      n_sampled = n_per, n_sero = round(p_inter * n_per)
    )
  }
  animal_periods <- do.call(rbind, ap_rows)

  d <- spillover_data(
    human = human, animal = animal, animal_indiv = animal_indiv,
    site = site_df, n_hosts = K, host_names = host_names
  )
  d$animal_periods <- animal_periods
  d
}

# ---- CASE 1: Sin Nombre hantavirus, USA (1993) -------------------------------
# Refs: Nichol et al. 1993 Science; Childs et al. 1994 J Infect Dis;
# Mills et al. 1997 J Mammal; Graham & Chomel 1997 EID.
# Published values: ~30.4% seroprev in P. maniculatus (n=813 across 21 sites);
# <5% in other rodents; endemic age accumulation; stable across years.
case_sin_nombre <- function() {
  n_sites <- 21
  host_names <- c("P. maniculatus", "P. truei", "Microtus", "Other rodents")
  prev_pcr <- list(
    "P. maniculatus" = pmax(0, rnorm(n_sites, 0.15, 0.08)),
    "P. truei"      = pmax(0, rnorm(n_sites, 0.03, 0.02)),
    "Microtus"      = pmax(0, rnorm(n_sites, 0.01, 0.01)),
    "Other rodents" = pmax(0, rnorm(n_sites, 0.005, 0.005))
  )
  prev_sero <- list(
    "P. maniculatus" = pmin(pmax(0, rnorm(n_sites, 0.304, 0.12)), 0.6),
    "P. truei"      = pmax(0, rnorm(n_sites, 0.04, 0.03)),
    "Microtus"      = pmax(0, rnorm(n_sites, 0.015, 0.01)),
    "Other rodents" = pmax(0, rnorm(n_sites, 0.008, 0.006))
  )
  n_sampled <- list(
    "P. maniculatus" = round(rnorm(n_sites, 40, 12)),
    "P. truei"      = round(rnorm(n_sites, 12, 4)),
    "Microtus"      = round(rnorm(n_sites, 8, 3)),
    "Other rodents" = round(rnorm(n_sites, 5, 2))
  )
  exposure_probs <- list(
    "P. maniculatus" = rep(0.55, n_sites), "P. truei" = rep(0.10, n_sites),
    "Microtus" = rep(0.08, n_sites), "Other rodents" = rep(0.05, n_sites)
  )
  exposure_rr <- list(
    "P. maniculatus" = rep(8.0, n_sites), "P. truei" = rep(1.1, n_sites),
    "Microtus" = rep(1.0, n_sites), "Other rodents" = rep(1.0, n_sites)
  )
  # Endemic rodent reservoir: stable inter-epidemic seropositivity
  interepi_ratio <- list(
    "P. maniculatus" = 0.85, "P. truei" = 0.20,
    "Microtus" = 0.15, "Other rodents" = 0.10
  )
  endemic_str <- list(
    "P. maniculatus" = 0.90, "P. truei" = 0.40,
    "Microtus" = 0.30, "Other rodents" = 0.25
  )
  build_outbreak_dataset(
    n_sites = n_sites, host_names = host_names,
    host_prev_pcr = prev_pcr, host_prev_sero = prev_sero,
    host_n_sampled = n_sampled, human_per_site = 8,
    human_exposure_probs = exposure_probs, spillover_intensity = 0.001,
    exposure_rr = exposure_rr,
    interepi_prev_ratio = interepi_ratio,
    endemic_signature_str = endemic_str
  )
}

# ---- CASE 2: Nipah virus, Malaysia (1998-99) ---------------------------------
# Refs: Chua et al. 2000 Science; Yob et al. 2001 EID; Goh et al. 2000 NEJM;
# Epstein et al. 2008 IJID; Pulliam et al. 2012 J R Soc Interface.
# Published values: 265 cases, 41% CFR, 93% pig contact among cases;
# pigs ~85% seroprev on affected farms during outbreak, near-zero
# inter-epidemic after culling; Pteropus 5-30% widespread, stable.
case_nipah <- function() {
  n_sites <- 15
  host_names <- c("Pigs", "Pteropus", "Equines", "Dogs")
  prev_pcr <- list(
    Pigs     = c(rep(0.70, 7), rep(0.05, 8)),
    Pteropus = pmax(0, rnorm(n_sites, 0.04, 0.02)),
    Equines  = pmax(0, rnorm(n_sites, 0.005, 0.005)),
    Dogs     = pmax(0, rnorm(n_sites, 0.002, 0.003))
  )
  prev_sero <- list(
    Pigs     = c(rep(0.85, 7), rep(0.05, 8)),  # outbreak vs unaffected
    Pteropus = c(rep(0.18, 5), rep(0.14, 5), rep(0.12, 5)),  # broad, stable
    Equines  = pmax(0, rnorm(n_sites, 0.02, 0.02)),
    Dogs     = pmax(0, rnorm(n_sites, 0.01, 0.01))
  )
  n_sampled <- list(
    Pigs     = round(rnorm(n_sites, 55, 15)),
    Pteropus = round(rnorm(n_sites, 40, 12)),
    Equines  = round(rnorm(n_sites, 12, 4)),
    Dogs     = round(rnorm(n_sites, 10, 3))
  )
  exposure_probs <- list(
    Pigs     = c(rep(0.85, 7), rep(0.10, 8)),
    Pteropus = rep(0.05, n_sites),
    Equines  = rep(0.10, n_sites),
    Dogs     = rep(0.40, n_sites)
  )
  exposure_rr <- list(
    Pigs     = rep(10.0, n_sites), Pteropus = rep(1.2, n_sites),
    Equines  = rep(1.0, n_sites),  Dogs     = rep(1.0, n_sites)
  )
  interepi_ratio <- list(
    Pigs = 0.03, Pteropus = 0.92, Equines = 0.30, Dogs = 0.25
  )
  endemic_str <- list(
    Pigs = 0.10, Pteropus = 0.92, Equines = 0.40, Dogs = 0.35
  )
  build_outbreak_dataset(
    n_sites = n_sites, host_names = host_names,
    host_prev_pcr = prev_pcr, host_prev_sero = prev_sero,
    host_n_sampled = n_sampled, human_per_site = 15,
    human_exposure_probs = exposure_probs, spillover_intensity = 0.002,
    exposure_rr = exposure_rr,
    interepi_prev_ratio = interepi_ratio,
    endemic_signature_str = endemic_str
  )
}

# ---- CASE 3: Rift Valley fever, Kenya (2006-07) ------------------------------
# Refs: Anyangu et al. 2010 AJTMH; LaBeaud et al. 2011 EID;
# Britch et al. 2013 PLoS ONE; Murithi et al. 2011 EpiInfect.
# Published values: 392 cases, 26% CFR; wild ungulates 31.8% sero (2007) vs
# 3.3% pre-outbreak (ratio 0.10); cattle 60% during outbreak, near-zero after.
case_rvf <- function() {
  n_sites <- 12
  host_names <- c("Cattle", "Wild ungulates", "Camels", "Small ruminants")
  prev_pcr <- list(
    "Cattle"          = c(rep(0.20, 6), rep(0.02, 6)),
    "Wild ungulates"  = pmax(0, rnorm(n_sites, 0.05, 0.03)),
    "Camels"          = c(rep(0.10, 3), rep(0.01, 9)),
    "Small ruminants" = pmax(0, rnorm(n_sites, 0.03, 0.02))
  )
  prev_sero <- list(
    "Cattle"          = c(rep(0.60, 6), rep(0.08, 6)),
    "Wild ungulates"  = pmin(pmax(0, rnorm(n_sites, 0.318, 0.10)), 0.55),
    "Camels"          = c(rep(0.35, 3), rep(0.04, 9)),
    "Small ruminants" = pmax(0, rnorm(n_sites, 0.08, 0.05))
  )
  n_sampled <- list(
    "Cattle"          = round(rnorm(n_sites, 60, 15)),
    "Wild ungulates"  = round(rnorm(n_sites, 35, 10)),
    "Camels"          = round(rnorm(n_sites, 18, 6)),
    "Small ruminants" = round(rnorm(n_sites, 40, 10))
  )
  exposure_probs <- list(
    "Cattle"          = rep(0.55, n_sites),
    "Wild ungulates"  = rep(0.06, n_sites),
    "Camels"          = c(rep(0.20, 3), rep(0.02, 9)),
    "Small ruminants" = rep(0.30, n_sites)
  )
  exposure_rr <- list(
    "Cattle"          = rep(7.0, n_sites),
    "Wild ungulates"  = rep(1.5, n_sites),
    "Camels"          = rep(1.5, n_sites),
    "Small ruminants" = rep(2.0, n_sites)
  )
  # Wild ungulates persist inter-epidemic (vertical mosquito transmission)
  # Cattle collapse (vaccination, short Ab decay)
  interepi_ratio <- list(
    "Cattle" = 0.08, "Wild ungulates" = 0.85,
    "Camels" = 0.40, "Small ruminants" = 0.20
  )
  endemic_str <- list(
    "Cattle" = 0.15, "Wild ungulates" = 0.85,
    "Camels" = 0.55, "Small ruminants" = 0.40
  )
  build_outbreak_dataset(
    n_sites = n_sites, host_names = host_names,
    host_prev_pcr = prev_pcr, host_prev_sero = prev_sero,
    host_n_sampled = n_sampled, human_per_site = 15,
    human_exposure_probs = exposure_probs, spillover_intensity = 0.002,
    exposure_rr = exposure_rr,
    interepi_prev_ratio = interepi_ratio,
    endemic_signature_str = endemic_str
  )
}


# ---- Run all three cases -----------------------------------------------------
run_case <- function(name, case_data, true_host) {
  cat("\n----", name, "----\n")
  cat("True reservoir (molecular):", true_host, "\n")
  res <- compute_pi(case_data, bootstrap = 100, verbose = FALSE)
  pi_vec <- res$pi; rank_vec <- res$rank
  df <- data.frame(host = names(pi_vec), Pi = round(pi_vec, 3), rank = rank_vec)
  df <- df[order(df$rank), ]
  print(df, row.names = FALSE)
  if (!is.null(res$bootstrap)) {
    ci <- t(apply(res$bootstrap, 2, stats::quantile,
                  probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
    cat("Bootstrap 95% CI:\n"); print(round(ci, 3))
  }
  list(case = name, data = case_data, result = res, true_host = true_host)
}

cat("== Building datasets ==\n")
data_sn <- case_sin_nombre()
data_nv <- case_nipah()
data_rv <- case_rvf()

res_sn <- run_case("Sin Nombre (USA, 1993)", data_sn, "P. maniculatus")
res_nv <- run_case("Nipah (Malaysia, 1998-99)", data_nv, "Pteropus")
res_rv <- run_case("Rift Valley fever (Kenya, 2006-07)", data_rv, "Wild ungulates")

# Save R results
saveRDS(list(sn = res_sn, nv = res_nv, rv = res_rv),
        file = "../data/case_results.rds")

# Export to JSON for Python figure pipeline
to_json_block <- function(r) {
  mrg <- margin_boot(r)
  list(
    case = r$case,
    true_host = r$true_host,
    hosts = names(r$result$pi),
    pi = as.numeric(r$result$pi),
    rank = as.integer(r$result$rank),
    criteria = as.data.frame(r$result$criteria),
    weights = as.numeric(r$result$weights_used),
    weight_names = names(r$result$weights_used),
    ci_low = as.numeric(apply(r$result$bootstrap, 2, stats::quantile, 0.025, na.rm = TRUE)),
    ci_med = as.numeric(apply(r$result$bootstrap, 2, stats::quantile, 0.500, na.rm = TRUE)),
    ci_high = as.numeric(apply(r$result$bootstrap, 2, stats::quantile, 0.975, na.rm = TRUE)),
    margin = mrg$point,
    margin_ci_low = mrg$lo,
    margin_ci_high = mrg$hi,
    margin_competitor = mrg$competitor
  )
}

# Discrimination margin: Pi of the true host minus the highest Pi among the
# remaining candidates. Its interval is the percentile interval of that
# DIFFERENCE across the same bootstrap replicates (Methods), which is not the
# difference of the two marginal percentile intervals: the latter ignores the
# correlation between hosts within a replicate and is far too wide.
margin_boot <- function(r) {
  pi_vec <- r$result$pi
  hosts  <- names(pi_vec)
  ti     <- which(hosts == r$true_host)
  oi     <- setdiff(seq_along(hosts), ti)
  bm     <- r$result$bootstrap
  d      <- bm[, ti] - apply(bm[, oi, drop = FALSE], 1, max)
  q      <- stats::quantile(d, c(0.025, 0.975), na.rm = TRUE)
  list(point      = as.numeric(pi_vec[ti] - max(pi_vec[oi])),
       lo         = as.numeric(q[1]),
       hi         = as.numeric(q[2]),
       competitor = hosts[oi][which.max(pi_vec[oi])])
}

out <- list(
  sn = to_json_block(res_sn),
  nv = to_json_block(res_nv),
  rv = to_json_block(res_rv)
)
write_json(out, "../data/case_results.json", auto_unbox = TRUE, pretty = TRUE, digits = 4)

cat("\n== DONE ==\n")
cat("Results: case_results.rds and case_results.json\n")
