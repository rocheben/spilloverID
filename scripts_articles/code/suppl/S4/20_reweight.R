# ============================================================================
# 20_reweight.R -- evidence-quality reweighting of Pi, C9-as-NA variant,
# and tipping-point sweeps. Consumes work/reliability_inputs.rds.
#
# Reliability model:  r_c = a_c * e_c * u_c * s_c   in [0,1]
#   a_c  availability / discriminative power (verifiable: all-NA? zero variance?)
#   e_c  evidence-class reliability = sum_class phi_{c,class} * rho_class
#        with phi from the mechanical field audit of each criterion's source
#   u_c  support saturation = n_clust_c / (n_clust_c + n0)
#        (n_clust = number of INDEPENDENT sampling clusters feeding criterion c)
#   s_c  serological-specificity robustness = 1 / (1 + kappa * delta_c)
#        delta_c = MEASURED mean |criterion displacement| at specificity 0.95
# Applied weights: w'_c propto w_c * r_c, renormalized over available criteria.
# ============================================================================
RI <- readRDS("reliability_inputs.rds")
source("../_load_cases.R")
CASES <- lapply(names(CASE_LAB), function(k) list(label = unname(CASE_LAB[k]), data = NULL, true_host = unname(CASE_TRUE[k])))
names(CASES) <- names(CASE_LAB)
.blt <- build_cases()
for (k in names(CASES)) CASES[[k]]$data <- .blt[[k]]
BASE <- lapply(CASES, function(cs) suppressWarnings(compute_pi(cs$data, bootstrap = 0)))
CRIT <- lapply(BASE, function(b) b$criteria)

# Basis = the criterion matrices of the three reconstructions,
# i.e. the criterion matrices the manuscript reports, rebuilt here by
# run_cases.R from the same constructors and the same seed.

share <- RI$share; sens <- RI$sens; SUP <- RI$support
W0 <- default_weights()
CRITNM <- names(W0)
CASENM <- c("sn", "nv", "rv")
TRUE_HOST <- sapply(CASES, `[[`, "true_host")

PROFILES <- list(
  conservative = list(rho = c(pcr = 1.00, sero = 0.45, decl = 0.35),
                      kappa = 3.0, n0 = 8),
  intermediate = list(rho = c(pcr = 1.00, sero = 0.70, decl = 0.60),
                      kappa = 1.5, n0 = 4)
)

# ---- availability / discriminative term, per case --------------------------
avail_term <- function(cm) {
  sapply(CRITNM, function(cn) {
    x <- cm[, cn]
    if (all(is.na(x))) return(0)              # criterion returns nothing
    if (stats::sd(x, na.rm = TRUE) == 0) return(0)  # constant => no ranking info
    1
  })
}

build_r <- function(nm, prof) {
  cm <- RI$CRIT[[nm]]
  sup <- SUP[[nm]]
  a <- avail_term(cm)
  e <- as.numeric(share[CRITNM, "phi_pcr"] * prof$rho["pcr"] +
                  share[CRITNM, "phi_sero"] * prof$rho["sero"] +
                  share[CRITNM, "phi_decl"] * prof$rho["decl"])
  names(e) <- CRITNM
  nc <- stats::setNames(sup$n_clust[match(CRITNM, sup$crit)], CRITNM)
  u <- nc / (nc + prof$n0)
  dl <- sens[nm, CRITNM]; dl[is.na(dl)] <- 0
  s <- 1 / (1 + prof$kappa * dl)
  r <- a * e * u * s
  data.frame(crit = CRITNM, a = a, e = e, n_clust = nc, u = u,
             delta_spec = as.numeric(dl), s = s, r = as.numeric(r),
             row.names = NULL)
}

# ---- Pi under an arbitrary weight vector, from a FIXED criterion matrix -----
pi_from <- function(cm, w) {
  keep <- !apply(cm, 2, function(x) all(is.na(x)))
  w <- w[CRITNM]; w[!keep] <- 0
  if (sum(w) <= 0) return(stats::setNames(rep(NA_real_, nrow(cm)), rownames(cm)))
  w <- w / sum(w)
  cs <- cm; cs[is.na(cs)] <- 0
  stats::setNames(as.numeric(cs %*% w), rownames(cm))
}

margin_rank <- function(pv, truth) {
  ord <- sort(pv, decreasing = TRUE)
  comp <- names(ord)[names(ord) != truth][1]
  list(pi_true = unname(pv[truth]), competitor = comp,
       pi_comp = unname(pv[comp]), margin = unname(pv[truth] - pv[comp]),
       rank_true = unname(rank(-pv, ties.method = "min")[truth]))
}

# ============================ (2) REWEIGHTED Pi ==============================
rows <- list(); rtabs <- list()
for (nm in CASENM) {
  cm <- RI$CRIT[[nm]]; truth <- TRUE_HOST[[nm]]
  variants <- list(published = W0)
  for (pn in names(PROFILES)) {
    rt <- build_r(nm, PROFILES[[pn]]); rtabs[[paste(nm, pn)]] <- cbind(case = nm, profile = pn, rt)
    variants[[pn]] <- stats::setNames(W0 * rt$r, CRITNM)
  }
  # C9-as-NA variant (weights as reported, C9 dropped rather than dead-weighted)
  cm9 <- cm; cm9[, "C9"] <- NA_real_
  for (vn in names(variants)) {
    mr <- margin_rank(pi_from(cm, variants[[vn]]), truth)
    rows[[length(rows) + 1]] <- data.frame(
      case = CASES[[nm]]$label, true_host = truth, weighting = vn,
      pi_true = mr$pi_true, competitor = mr$competitor, pi_comp = mr$pi_comp,
      margin = mr$margin, rank_true = mr$rank_true, C9_treatment = "as reported (0)")
  }
  mr <- margin_rank(pi_from(cm9, W0), truth)
  rows[[length(rows) + 1]] <- data.frame(
    case = CASES[[nm]]$label, true_host = truth, weighting = "published",
    pi_true = mr$pi_true, competitor = mr$competitor, pi_comp = mr$pi_comp,
    margin = mr$margin, rank_true = mr$rank_true, C9_treatment = "set to NA")
  for (pn in names(PROFILES)) {
    mr <- margin_rank(pi_from(cm9, variants[[pn]]), truth)
    rows[[length(rows) + 1]] <- data.frame(
      case = CASES[[nm]]$label, true_host = truth, weighting = pn,
      pi_true = mr$pi_true, competitor = mr$competitor, pi_comp = mr$pi_comp,
      margin = mr$margin, rank_true = mr$rank_true, C9_treatment = "set to NA")
  }
}
TAB <- do.call(rbind, rows)
RTAB <- do.call(rbind, rtabs)

# effective weight actually carried by the maintenance triad, per variant
eff_triad <- list()
for (nm in CASENM) {
  cm <- RI$CRIT[[nm]]
  keep <- !apply(cm, 2, function(x) all(is.na(x)))
  for (vn in c("published", names(PROFILES))) {
    w <- if (vn == "published") W0 else W0 * build_r(nm, PROFILES[[vn]])$r
    w[!keep] <- 0; w <- w / sum(w)
    eff_triad[[length(eff_triad) + 1]] <- data.frame(
      case = nm, weighting = vn,
      w_C12 = w["C12"], w_C13 = w["C13"], w_C14 = w["C14"],
      w_triad = sum(w[c("C12", "C13", "C14")]),
      w_sero_weighted = sum(w * share[CRITNM, "phi_sero"]))
  }
}
EFF <- do.call(rbind, eff_triad); rownames(EFF) <- NULL

# ============================ (3) TIPPING POINTS =============================
QGRID <- seq(0, 1, by = 0.002)
disc_weights <- function(w, target, q) {
  if (identical(target, "sero")) {
    # discount proportional to each criterion's serological content
    w * (1 - (1 - q) * share[CRITNM, "phi_sero"])
  } else {
    w[target] <- w[target] * q; w
  }
}
TARGETS <- list(triad = c("C12", "C13", "C14"), C12 = "C12", sero = "sero")

sweep_rows <- list()
for (nm in CASENM) for (vn in c("published", names(PROFILES))) {
  cm <- RI$CRIT[[nm]]; truth <- TRUE_HOST[[nm]]
  wbase <- if (vn == "published") W0 else stats::setNames(W0 * build_r(nm, PROFILES[[vn]])$r, CRITNM)
  for (tn in names(TARGETS)) {
    m <- sapply(QGRID, function(q)
      margin_rank(pi_from(cm, disc_weights(wbase, TARGETS[[tn]], q)), truth)$margin)
    sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
      case = nm, case_label = CASES[[nm]]$label, weighting = vn,
      target = tn, q = QGRID, margin = m)
  }
}
SWEEP <- do.call(rbind, sweep_rows)

# tipping point q*: largest q at which margin first becomes <= 0 scanning down
tip <- function(df) {
  df <- df[order(-df$q), ]
  i <- which(df$margin <= 0)[1]
  if (is.na(i)) return(NA_real_)
  if (i == 1) return(df$q[1])
  q_hi <- df$q[i - 1]; q_lo <- df$q[i]
  m_hi <- df$margin[i - 1]; m_lo <- df$margin[i]
  q_lo + (0 - m_lo) * (q_hi - q_lo) / (m_hi - m_lo)   # linear interpolation
}
TIP <- do.call(rbind, lapply(split(SWEEP, list(SWEEP$case, SWEEP$weighting, SWEEP$target), drop = TRUE),
  function(df) data.frame(case = df$case[1], case_label = df$case_label[1],
                          weighting = df$weighting[1], target = df$target[1],
                          margin_q1 = df$margin[df$q == 1], margin_q0 = df$margin[df$q == 0],
                          q_star = tip(df))))
rownames(TIP) <- NULL
# Interpretation guard: q = 1 is the undiscounted baseline. If the margin is
# already negative there, the rank is flipped BEFORE any discount is applied and
# q* is not a tipping point -- flag it rather than reporting a spurious 1.00.
TIP$status <- ifelse(TIP$margin_q1 <= 0, "already mis-ranked at baseline",
               ifelse(is.na(TIP$q_star), "never flips (holds down to q=0)", "flips at q*"))
TIP$q_star[TIP$margin_q1 <= 0] <- NA_real_

dir.create("out", showWarnings = FALSE)
write.csv(TAB,   "../../../data/suppl_S4_reliability_weighting/table_reweighted_pi.csv", row.names = FALSE)
write.csv(RTAB,  "../../../data/suppl_S4_reliability_weighting/table_reliability_factors.csv", row.names = FALSE)
write.csv(EFF,   "../../../data/suppl_S4_reliability_weighting/table_effective_weights.csv", row.names = FALSE)
write.csv(TIP,   "../../../data/suppl_S4_reliability_weighting/table_tipping_points.csv", row.names = FALSE)
write.csv(SWEEP, "../../../data/suppl_S4_reliability_weighting/sweep_margin_vs_q.csv", row.names = FALSE)

cat("\n=== reliability factors r_c (conservative, Nipah) ===\n")
print(subset(RTAB, case == "nv" & profile == "conservative")[, c("crit","a","e","n_clust","u","delta_spec","s","r")],
      row.names = FALSE, digits = 3)
cat("\n=== effective weights ===\n"); print(EFF, row.names = FALSE, digits = 3)
cat("\n=== reweighted Pi table ===\n"); print(TAB, row.names = FALSE, digits = 3)
cat("\n=== tipping points ===\n")
print(TIP[order(TIP$target, TIP$case, TIP$weighting), c("case","weighting","target","margin_q1","margin_q0","q_star")],
      row.names = FALSE, digits = 3)
