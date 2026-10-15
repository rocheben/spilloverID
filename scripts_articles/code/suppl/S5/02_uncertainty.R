# ============================================================================
# Uncertainty propagation for the composite plausibility score Pi.
#
# The published bootstrap (compute_pi(..., bootstrap = B)) resamples SITES with
# replacement via .resample_data(), which carries data$human, data$animal,
# data$animal_indiv and data$site along. It never touches data$animal_periods,
# the sole support of C12 -- the highest-weighted criterion (w = 0.174 effective
# in the three real cases). In all three case datasets animal_periods contains
# a single pseudo-site (site = 1), so even a site-aware resampler would leave
# C12 exactly constant across replicates.
#
# This script re-runs the three cases under four resampling schemes:
#   published : .resample_data() verbatim                        (site cluster)
#   c12       : site cluster + binomial resample of animal_periods
#   c13       : site cluster + within-site resample of animal_indiv
#   corrected : site cluster + both
# and reports interval widths and P(rank 1) per host.
# ============================================================================
source("../_load_cases.R")
dir.create("../../../data/suppl_S5_uncertainty", showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(parallel))

B      <- as.integer(Sys.getenv("S5_B", "300"))
NCORES <- { n <- suppressWarnings(detectCores()); if (!is.finite(n)) n <- as.integer(Sys.getenv("NPROC", "2")); max(1L, as.integer(Sys.getenv("S5_N_CORE", as.character(n)))) }

# ---- resamplers -------------------------------------------------------------

# Verbatim copy of the published .resample_data (site-level cluster bootstrap),
# plus optional additions. `indiv` = within-site resampling of animal_indiv
# rows; `periods` = binomial resampling of animal_periods seropositive counts.
resample_data2 <- function(data, sampled_sites, indiv = FALSE, periods = FALSE) {
  new_human <- list(); new_animal <- list(); new_animal_indiv <- list(); new_site <- list()
  for (i in seq_along(sampled_sites)) {
    s <- sampled_sites[i]; new_id <- i
    h <- data$human[data$human$site == s, ]; h$site <- new_id; new_human[[i]] <- h
    a <- data$animal[data$animal$site == s, ]; a$site <- new_id; new_animal[[i]] <- a
    if (!is.null(data$animal_indiv)) {
      ai <- data$animal_indiv[data$animal_indiv$site == s, ]
      if (indiv && nrow(ai) > 0) {
        # stratify the within-site resample by host so that per-host sample
        # sizes (which drive the C13 catalytic fit) are preserved
        ai <- do.call(rbind, lapply(split(ai, ai$host), function(g)
          g[sample.int(nrow(g), nrow(g), replace = TRUE), , drop = FALSE]))
      }
      ai$site <- new_id
      new_animal_indiv[[i]] <- ai
    }
    sr <- data$site[data$site$site == s, ]; sr$site <- new_id; new_site[[i]] <- sr
  }
  data2 <- data
  data2$human <- do.call(rbind, new_human)
  data2$animal <- do.call(rbind, new_animal)
  if (!is.null(data$animal_indiv)) data2$animal_indiv <- do.call(rbind, new_animal_indiv)
  data2$site <- do.call(rbind, new_site)
  if (periods && !is.null(data$animal_periods)) {
    ap <- data$animal_periods
    p_hat <- pmin(pmax(ap$n_sero / pmax(ap$n_sampled, 1), 0), 1)
    ap$n_sero <- rbinom(nrow(ap), size = ap$n_sampled, prob = p_hat)
    data2$animal_periods <- ap
  }
  data2
}

MODES <- list(
  published = list(indiv = FALSE, periods = FALSE),
  c12       = list(indiv = FALSE, periods = TRUE),
  c13       = list(indiv = TRUE,  periods = FALSE),
  corrected = list(indiv = TRUE,  periods = TRUE)
)

# ---- one bootstrap replicate ------------------------------------------------
boot_rep <- function(b, data, mode, seed_base) {
  set.seed(seed_base + b)
  sites <- unique(data$human$site)
  sampled <- sample(sites, replace = TRUE)
  d_b <- resample_data2(data, sampled, indiv = mode$indiv, periods = mode$periods)
  r <- tryCatch(suppressWarnings(compute_pi(d_b, bootstrap = 0)),
                error = function(e) NULL)
  if (is.null(r)) return(rep(NA_real_, data$n_hosts))
  as.numeric(r$pi)
}

run_mode <- function(data, mode, seed_base) {
  reps <- mclapply(seq_len(B), boot_rep, data = data, mode = mode,
                   seed_base = seed_base, mc.cores = NCORES)
  m <- do.call(rbind, reps)
  colnames(m) <- data$host_names
  m
}

# ---- main -------------------------------------------------------------------
dat <- build_cases()
point <- lapply(dat, function(d) suppressWarnings(compute_pi(d, bootstrap = 0)))

rows <- list(); boots <- list()
for (nm in names(dat)) {
  d <- dat[[nm]]
  cat("== case", nm, "==\n"); flush.console()
  for (mo in names(MODES)) {
    t0 <- Sys.time()
    bm <- run_mode(d, MODES[[mo]], seed_base = 771000 + 1000 * match(nm, names(dat)) +
                     100 * match(mo, names(MODES)))
    boots[[paste(nm, mo, sep = "|")]] <- bm
    ok <- rowSums(is.na(bm)) == 0
    bmo <- bm[ok, , drop = FALSE]
    # P(rank 1): fraction of replicates in which the host has the strictly
    # highest Pi (ties split evenly among the tied hosts)
    rk1 <- t(apply(bmo, 1, function(x) { mx <- max(x); as.numeric(x == mx) / sum(x == mx) }))
    for (k in seq_len(d$n_hosts)) {
      q <- quantile(bmo[, k], c(0.025, 0.25, 0.5, 0.75, 0.975))
      rows[[length(rows) + 1]] <- data.frame(
        case = nm, case_label = unname(CASE_LAB[nm]), host = d$host_names[k],
        is_true_host = d$host_names[k] == unname(CASE_TRUE[nm]),
        scheme = mo, n_valid = nrow(bmo),
        pi_point = as.numeric(point[[nm]]$pi[k]),
        ci_low = as.numeric(q[1]), ci_q25 = as.numeric(q[2]),
        ci_med = as.numeric(q[3]), ci_q75 = as.numeric(q[4]),
        ci_high = as.numeric(q[5]),
        width = as.numeric(q[5] - q[1]),
        boot_sd = sd(bmo[, k]),
        p_rank1 = mean(rk1[, k]),
        stringsAsFactors = FALSE)
    }
    cat(sprintf("   %-10s %5.1fs  mean width = %.3f\n", mo,
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                mean(sapply(seq_len(d$n_hosts), function(k)
                  diff(quantile(bmo[, k], c(0.025, 0.975))))))); flush.console()
  }
}
res <- do.call(rbind, rows)

# ---- margin (top1 - top2) distribution per case and scheme ------------------
mrows <- list()
for (key in names(boots)) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  nm <- parts[1]; mo <- parts[2]
  bm <- boots[[key]]; bm <- bm[rowSums(is.na(bm)) == 0, , drop = FALSE]
  srt <- t(apply(bm, 1, sort, decreasing = TRUE))
  marg <- srt[, 1] - srt[, 2]
  ihat <- match(unname(CASE_TRUE[nm]), colnames(bm))
  marg_true <- bm[, ihat] - apply(bm[, -ihat, drop = FALSE], 1, max)
  q <- quantile(marg_true, c(0.025, 0.5, 0.975))
  mrows[[length(mrows) + 1]] <- data.frame(
    case = nm, case_label = unname(CASE_LAB[nm]), scheme = mo,
    margin_point = as.numeric(sort(point[[nm]]$pi, decreasing = TRUE)[1] -
                              sort(point[[nm]]$pi, decreasing = TRUE)[2]),
    margin_true_med = as.numeric(q[2]), margin_true_lo = as.numeric(q[1]),
    margin_true_hi = as.numeric(q[3]),
    margin_true_se = sd(marg_true),
    p_true_rank1 = mean(marg_true > 0) + 0.5 * mean(marg_true == 0),
    stringsAsFactors = FALSE)
}
margins <- do.call(rbind, mrows)

# ---- diagnostic: is C12 constant under the published bootstrap? -------------
diag_rows <- list()
for (nm in names(dat)) {
  d <- dat[[nm]]
  for (mo in c("published", "corrected")) {
    cs <- mclapply(1:100, function(b) {
      set.seed(99000 + b)
      sampled <- sample(unique(d$human$site), replace = TRUE)
      d_b <- resample_data2(d, sampled, indiv = MODES[[mo]]$indiv,
                            periods = MODES[[mo]]$periods)
      r <- tryCatch(suppressWarnings(compute_pi(d_b, bootstrap = 0)), error = function(e) NULL)
      if (is.null(r)) return(rep(NA_real_, 2 * d$n_hosts))
      c(r$criteria[, "C12"], r$criteria[, "C13"])
    }, mc.cores = NCORES)
    M <- do.call(rbind, cs); M <- M[rowSums(is.na(M)) == 0, , drop = FALSE]
    K <- d$n_hosts
    diag_rows[[length(diag_rows) + 1]] <- data.frame(
      case = nm, scheme = mo,
      C12_sd_mean = mean(apply(M[, 1:K, drop = FALSE], 2, sd)),
      C12_n_distinct = mean(apply(M[, 1:K, drop = FALSE], 2, function(x) length(unique(round(x, 8))))),
      C13_sd_mean = mean(apply(M[, (K+1):(2*K), drop = FALSE], 2, sd)),
      C13_n_distinct = mean(apply(M[, (K+1):(2*K), drop = FALSE], 2, function(x) length(unique(round(x, 8))))),
      stringsAsFactors = FALSE)
  }
}
crit_diag <- do.call(rbind, diag_rows)

dir.create("out", showWarnings = FALSE)
write.csv(res,       "../../../data/suppl_S5_uncertainty/pi_intervals_by_case_host.csv", row.names = FALSE)
write.csv(margins,   "../../../data/suppl_S5_uncertainty/pi_margins_by_case.csv",        row.names = FALSE)
write.csv(crit_diag, "../../../data/suppl_S5_uncertainty/c12_c13_bootstrap_variability.csv", row.names = FALSE)
saveRDS(list(boots = boots, point = point, B = B), "../../../data/suppl_S5_uncertainty/boot_raw.rds")

cat("\n--- interval width, published vs corrected ---\n")
w <- reshape(res[, c("case","host","scheme","width")], idvar = c("case","host"),
             timevar = "scheme", direction = "wide")
w$delta <- w$width.corrected - w$width.published
print(w, row.names = FALSE, digits = 3)
cat("\n--- margin / P(rank 1) ---\n")
print(margins, row.names = FALSE, digits = 3)
cat("\n--- C12 / C13 variability across bootstrap replicates ---\n")
print(crit_diag, row.names = FALSE, digits = 3)
cat("\nDONE 02\n")
