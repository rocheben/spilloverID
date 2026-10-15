#' ============================================================================
#' x10_baselines.R -- External comparators / reference rankings for Pi
#'
#' Reference rankings against which the composite Pi is scored, all computed
#' from exactly the same spillover_data object that Pi reads.
#'
#' An obvious comparator would be the Mollentze & Streicker (2020 PNAS) /
#' Becker et al. (2022 Lancet Microbe) family of zoonotic-risk predictors.
#' They cannot be run here, and we do not fabricate scores from them. They take as input a viral genome (dinucleotide/CpG composition,
#' codon usage) or host-level phylogenetic and life-history traits, and return a
#' PRIOR probability that a virus species infects humans, or that a mammal
#' species is a reservoir *somewhere*. They are fitted across virus/host species
#' pools BEFORE any outbreak, they do not read outbreak observations, and every
#' host in a within-outbreak candidate set already has a fixed species-level
#' score that is identical for every outbreak that species is ever involved in.
#' Running them here would require host genomes and the published model objects,
#' neither of which is in this bundle, and their output would not vary with the
#' data Pi consumes. We therefore do NOT fabricate scores from those papers.
#' Instead B6 below is an explicitly labelled CARICATURE of the *only*
#' transferable part of that idea -- an a-priori taxonomic-order prior.
#'
#' Baselines implemented (all read the SAME spillover_data that Pi reads):
#'   B1 sero      pooled animal seroprevalence per host
#'   B2 pcr       pooled animal PCR prevalence per host
#'   B3 or        crude exposure odds ratio (declared exposure vs human sero)
#'   B4 C12       inter-epidemic persistence criterion, alone
#'   B5 C7        best single criterion on the simulated benchmark
#'                (see fig_compindiv_sim.csv)
#'   B6 taxo      a-priori taxonomic-order prior  [CARICATURE PROXY]
#' ============================================================================

BASELINE_IDS <- c("pi", "sero", "pcr", "or", "C12", "C7", "C14", "C5", "taxo")

BASELINE_LABELS <- c(
  pi   = "Composite Pi (14 criteria)",
  sero = "Animal seroprevalence alone",
  pcr  = "Animal PCR prevalence alone",
  or   = "Exposure odds ratio (crude)",
  C12  = "C12 inter-epidemic persistence alone",
  C7   = "C7 catalytic FOI alone (best single, suppl.)",
  C14  = "C14 spatial independence alone",
  C5   = "C5 declared-exposure fraction alone",
  taxo = "A-priori taxonomic order [proxy caricature]"
)

# ---- B1 / B2: pooled animal prevalence --------------------------------------
# Sample-size weighted (sum n_sero / sum n_sampled), i.e. what a field team
# reports as "seroprevalence in species k across the survey".
baseline_prevalence <- function(data, what = c("sero", "pcr")) {
  what <- match.arg(what)
  K <- data$n_hosts
  a <- data$animal
  num_col <- if (what == "sero") "n_sero" else "n_pcr"
  out <- vapply(seq_len(K), function(k) {
    sub <- a[a$host == k, , drop = FALSE]
    if (nrow(sub) == 0) return(NA_real_)
    if (num_col %in% names(sub)) {
      den <- sum(sub$n_sampled, na.rm = TRUE)
      if (is.finite(den) && den > 0) return(sum(sub[[num_col]], na.rm = TRUE) / den)
    }
    mean(sub[[paste0("prev_", what)]], na.rm = TRUE)
  }, numeric(1))
  stats::setNames(out, data$host_names)
}

# Unweighted site mean, kept as a sensitivity variant.
baseline_prevalence_unweighted <- function(data, what = c("sero", "pcr")) {
  what <- match.arg(what)
  K <- data$n_hosts
  a <- data$animal
  out <- vapply(seq_len(K), function(k)
    mean(a[[paste0("prev_", what)]][a$host == k], na.rm = TRUE), numeric(1))
  stats::setNames(out, data$host_names)
}

# ---- B3: crude exposure odds ratio ------------------------------------------
# 2x2 table of declared exposure to host k vs human seropositivity, pooled over
# sites, Haldane-Anscombe (+0.5) continuity correction. Deliberately CRUDE: the
# site-adjusted mixed-model version of this quantity is already criterion C3.
baseline_exposure_or <- function(data) {
  K <- data$n_hosts
  d <- data$human
  out <- vapply(seq_len(K), function(k) {
    ec <- paste0("exp_h", k)
    if (!ec %in% names(d)) return(NA_real_)
    e <- d[[ec]]; y <- d$sero
    ok <- !is.na(e) & !is.na(y)
    e <- e[ok]; y <- y[ok]
    a <- sum(e == 1 & y == 1); b <- sum(e == 1 & y == 0)
    c_ <- sum(e == 0 & y == 1); dd <- sum(e == 0 & y == 0)
    ((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (c_ + 0.5))
  }, numeric(1))
  stats::setNames(out, data$host_names)
}

# ---- B6: a-priori taxonomic-order prior  [CARICATURE] ------------------------
# NOT a Mollentze & Streicker score. A deliberately crude stand-in for the one
# element of trait-based predictors that transfers to a within-outbreak setting:
# the fact that some mammalian orders are cited as reservoirs far more often
# than others. Tiers fixed a priori, before looking at any Pi output:
#   3 Chiroptera            2 Rodentia
#   1 wild non-volant non-rodent mammals      0 domestic animals
# Structural ties are intended and are the point: the proxy cannot separate four
# rodent species from one another. Ranking below therefore also reports the
# expected top-1 under random tie-breaking.
TAXO_TIERS <- c(
  # Chiroptera
  "Pteropus" = 3,
  # Rodentia
  "P. maniculatus" = 2, "P. truei" = 2, "Microtus" = 2, "Other rodents" = 2,
  # wild non-volant non-rodent mammals
  "Wild ungulates" = 1,
  # domestic
  "Pigs" = 0, "Cattle" = 0, "Camels" = 0, "Small ruminants" = 0,
  "Equines" = 0, "Dogs" = 0
)

baseline_taxonomic <- function(data) {
  nm <- data$host_names
  sc <- TAXO_TIERS[nm]
  # Unnamed / synthetic hosts (h1..h4 in the simulator) have no taxonomy at all:
  # the proxy is undefined, which is itself the finding.
  if (anyNA(sc)) return(stats::setNames(rep(NA_real_, length(nm)), nm))
  stats::setNames(as.numeric(sc), nm)
}

# ---- ranking helpers ---------------------------------------------------------
#' Rank a score vector, higher = better, and score the true reservoir.
#'
#' top1_strict follows the manuscript's own convention (comp_vs_indiv.R): the
#' true host must hold a UNIQUE maximum. top1_exp is the expected top-1 under
#' random tie-breaking (1/t when the true host is one of t tied maxima), which
#' is the fair convention for the deliberately tie-heavy taxonomic proxy.
rank_true_host <- function(scores, k_true) {
  K <- length(scores)
  if (is.na(k_true))
    return(list(rank = NA_integer_, top1_strict = NA_integer_,
                top1_exp = NA_real_, top2 = NA_integer_, n_tied_top = NA_integer_))
  if (all(is.na(scores)))
    return(list(rank = NA_integer_, top1_strict = NA_integer_,
                top1_exp = NA_real_, top2 = NA_integer_, n_tied_top = NA_integer_))
  s <- as.numeric(scores); s[is.na(s)] <- -Inf
  rk <- rank(-s, ties.method = "min")
  n_tied_top <- sum(s == max(s))
  is_top_tied <- s[k_true] == max(s)
  list(
    rank        = as.integer(rk[k_true]),
    top1_strict = as.integer(rk[k_true] == 1 && n_tied_top == 1),
    top1_exp    = if (is_top_tied) 1 / n_tied_top else 0,
    top2        = as.integer(rk[k_true] <= 2),
    n_tied_top  = as.integer(n_tied_top)
  )
}

#' Scale-free confidence of the top pick, for the S11 negative control.
#' g = (s(1) - s(2)) / (s(1) - s(K)), in [0,1]; NA if the score is flat.
#' A method that produces a large g when there is nothing to find is
#' over-confident on the negative control.
top_gap <- function(scores) {
  s <- sort(as.numeric(scores[!is.na(scores)]), decreasing = TRUE)
  if (length(s) < 3) return(NA_real_)
  spread <- s[1] - s[length(s)]
  if (!is.finite(spread) || spread <= 0) return(NA_real_)
  (s[1] - s[2]) / spread
}

#' Compute every baseline score vector for one dataset.
#' `crit` is the criteria matrix from compute_pi() (hosts x C1..C14) so that
#' C5/C7/C12/C14 are byte-identical to the published criterion functions.
all_baseline_scores <- function(data, pi_vec, crit) {
  list(
    pi   = as.numeric(pi_vec),
    sero = as.numeric(baseline_prevalence(data, "sero")),
    pcr  = as.numeric(baseline_prevalence(data, "pcr")),
    or   = as.numeric(baseline_exposure_or(data)),
    C12  = as.numeric(crit[, "C12"]),
    C7   = as.numeric(crit[, "C7"]),
    C14  = as.numeric(crit[, "C14"]),
    C5   = as.numeric(crit[, "C5"]),
    taxo = as.numeric(baseline_taxonomic(data))
  )
}
