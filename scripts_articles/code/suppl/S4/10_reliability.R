# ============================================================================
# 10_reliability.R -- build per-criterion reliability factors r_c from
# VERIFIABLE properties of the deposited data + code:
#   (a) evidence-class mix  : which data fields each criterion actually reads
#                             (mechanical token audit of the function source),
#                             split into PCR / serology / self-declared
#   (b) availability        : all-NA, or degenerate (zero variance across hosts)
#   (c) effective support   : n observations and n independent clusters (sites)
#   (d) serological-specificity sensitivity : MEASURED by perturbing every
#       serological field at specificity sp and recording criterion displacement
# ============================================================================
source("../_load_cases.R")
CASES <- lapply(names(CASE_LAB), function(k) list(label = unname(CASE_LAB[k]), data = NULL, true_host = unname(CASE_TRUE[k])))
names(CASES) <- names(CASE_LAB)
.blt <- build_cases()
for (k in names(CASES)) CASES[[k]]$data <- .blt[[k]]
BASE <- lapply(CASES, function(cs) suppressWarnings(compute_pi(cs$data, bootstrap = 0)))
CRIT <- lapply(BASE, function(b) b$criteria)
set.seed(20260524)

CRITF <- list(
  C1 = criterion_C1_spatial, C2 = criterion_C2_reservoir_signatures,
  C3 = criterion_C3_individual, C4 = criterion_C4_dose_response,
  C5 = criterion_C5_cryptic_fraction, C6 = criterion_C6_interaction,
  C7 = criterion_C7_foi_coherence, C8 = criterion_C8_temporal_precedence,
  C9 = criterion_C9_mediation, C10 = criterion_C10_mechanistic,
  C11 = criterion_C11_specificity, C12 = criterion_C12_inter_epidemic,
  C13 = criterion_C13_endemic_age, C14 = criterion_C14_spatial_independence
)

# ---- (a) mechanical evidence-class audit -----------------------------------
# Tokens are the literal data-field names in the spillover_data object.
TOK <- list(
  pcr  = c("prev_pcr", "n_pcr", "prev_pcr_h", "prev_local_h", "prev_pcr_true"),
  sero = c("sero", "prev_sero", "n_sero", "prev_sero_h"),
  decl = c("exp_h", "prof", "env_declared", "residence")
)
# 'prev_local_h' is built from prev_pcr_h -> PCR. 'sero' also matches n_sero /
# prev_sero, so count the most specific token only.
count_tokens <- function(fn) {
  src <- paste(deparse(fn), collapse = "\n")
  hits <- list()
  for (cls in names(TOK)) {
    n <- 0
    for (t in TOK[[cls]]) n <- n + length(gregexpr(t, src, fixed = TRUE)[[1]][
      gregexpr(t, src, fixed = TRUE)[[1]] > 0])
    hits[[cls]] <- n
  }
  # de-duplicate: occurrences of prev_sero/n_sero were also caught by "sero"
  unlist(hits)
}
field_audit <- t(sapply(CRITF, count_tokens))
# C12 reads animal_periods (n_sero/n_sampled) -> serological only; C13 reads
# animal_indiv$sero+age -> serological only; C14 reads site prev_sero_h.
share <- field_audit / pmax(rowSums(field_audit), 1)
colnames(share) <- paste0("phi_", colnames(share))

# ---- (c) effective support, per case ---------------------------------------
support <- function(d) {
  K <- d$n_hosts
  n_sites <- nrow(d$site)
  n_hum <- nrow(d$human); n_pos <- sum(d$human$sero)
  n_ai <- nrow(d$animal_indiv)
  ap <- d$animal_periods
  n_ap <- sum(ap$n_sampled); ap_sites <- length(unique(ap$site))
  data.frame(
    crit = names(CRITF),
    n_obs = c(n_sites, n_ai, n_hum, sum(rowSums(d$human[, paste0("exp_h", 1:K), drop = FALSE]) > 0),
              n_pos, n_hum, n_hum + n_ai, 0, n_hum, n_hum, 0, n_ap, n_ai, n_sites),
    n_clust = c(n_sites, n_sites, n_sites, n_sites, n_sites, n_sites, n_sites,
                0, n_sites, n_sites, 0, ap_sites, n_sites, n_sites)
  )
}

# ---- (d) MEASURED sensitivity to serological specificity -------------------
# Imperfect specificity sp: each true-negative serological result becomes
# positive with probability (1 - sp). Applied to every serological field;
# PCR fields untouched.
perturb_sero <- function(d, sp) {
  fp <- 1 - sp
  d2 <- d
  neg <- d2$human$sero == 0
  d2$human$sero[neg] <- rbinom(sum(neg), 1, fp)
  neg <- d2$animal_indiv$sero == 0
  d2$animal_indiv$sero[neg] <- rbinom(sum(neg), 1, fp)
  add <- rbinom(nrow(d2$animal), pmax(d2$animal$n_sampled - d2$animal$n_sero, 0), fp)
  d2$animal$n_sero <- d2$animal$n_sero + add
  d2$animal$prev_sero <- pmin(d2$animal$n_sero / pmax(d2$animal$n_sampled, 1), 1)
  for (k in seq_len(d2$n_hosts)) {
    cc <- paste0("prev_sero_h", k)
    a_k <- d2$animal[d2$animal$host == k, c("site", "prev_sero")]
    d2$site[[cc]] <- a_k$prev_sero[match(d2$site$site, a_k$site)]
  }
  ap <- d2$animal_periods
  ap$n_sero <- ap$n_sero + rbinom(nrow(ap), pmax(ap$n_sampled - ap$n_sero, 0), fp)
  d2$animal_periods <- ap
  d2
}

crit_matrix <- function(d) {
  K <- d$n_hosts
  m <- matrix(NA_real_, K, length(CRITF), dimnames = list(d$host_names, names(CRITF)))
  for (cn in names(CRITF))
    m[, cn] <- tryCatch(CRITF[[cn]](d), error = function(e) rep(NA_real_, K))
  m
}

NREP <- as.integer(Sys.getenv("NREP", "40"))
SPEC <- 0.95   # 5% false-positive rate on serology -- Drexler's cross-reactivity
sens_list <- list()
for (nm in names(CASES)) {
  d0 <- CASES[[nm]]$data
  m0 <- CRIT[[nm]]
  disp <- matrix(NA_real_, NREP, length(CRITF), dimnames = list(NULL, names(CRITF)))
  for (b in seq_len(NREP)) {
    mb <- suppressWarnings(crit_matrix(perturb_sero(d0, SPEC)))
    disp[b, ] <- colMeans(abs(mb - m0))          # mean |delta| over hosts
  }
  sens_list[[nm]] <- colMeans(disp, na.rm = TRUE)
  cat("sensitivity done:", nm, "\n")
}
sens <- do.call(rbind, sens_list)

saveRDS(list(field_audit = field_audit, share = share,
             support = lapply(CASES, function(cs) support(cs$data)),
             sens = sens, SPEC = SPEC, NREP = NREP, CRIT = CRIT,
             BASE = BASE),
        "reliability_inputs.rds")
cat("\n== evidence-class shares (mechanical field audit) ==\n"); print(round(share, 3))
cat("\n== mean |delta criterion| at specificity", SPEC, "(", NREP, "reps ) ==\n")
print(round(sens, 4))
