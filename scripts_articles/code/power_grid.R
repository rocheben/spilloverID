#' ============================================================================
#' Prospective power analysis (Supplementary Figure 6)
#'
#' Top-1 identification rate as a function of the number of sampled sites and
#' the number of animals sampled per host per site, for the canonical scenario
#' S1 and the bridge-host scenario S2. Humans per site is held at the default
#' of 60 throughout, so the grid isolates the two design choices a field team
#' actually controls when planning animal sampling.
#'
#' Writes ../data/power_grid.json
#' Runtime: a few minutes.
#' ============================================================================

suppressPackageStartupMessages({
  library(spilloverID)
  library(jsonlite)
})

N_REP     <- 15
SITES     <- c(10, 15, 20, 30)
ANIMALS   <- c(10, 20, 40, 60)
SCENARIOS <- c("S1_canonical", "S2_bridge")
TARGET    <- 80  # top-1 target, in percent

rows <- list()
for (sc in SCENARIOS) {
  for (ns in SITES) {
    for (na in ANIMALS) {
      p <- scenario_params(sc)
      p$n_sites <- ns
      p$n_animals_per_host_site <- na
      hit <- integer(N_REP)
      for (rep in seq_len(N_REP)) {
        seed <- 20260524 + 17 * rep + 1000 * ns + 7 * na +
                100000 * match(sc, SCENARIOS)
        sim <- simulate_spillover(p, seed = seed)
        truth <- attr(sim, "truth")$k_star
        res <- suppressWarnings(compute_pi(sim, bootstrap = 0))
        hit[rep] <- as.integer(res$rank[truth] == 1)
      }
      rows[[length(rows) + 1]] <- data.frame(
        scenario = sc, n_sites = ns, n_animals = na,
        n_rep = N_REP, top1_pct = 100 * mean(hit),
        meets_target = 100 * mean(hit) >= TARGET
      )
      cat(sprintf("%-14s sites=%2d animals=%2d  top1=%5.1f%%\n",
                  sc, ns, na, 100 * mean(hit)))
    }
  }
}

grid <- do.call(rbind, rows)
write_json(list(target_pct = TARGET, n_rep_per_cell = N_REP,
                n_humans_per_site = scenario_params("S1_canonical")$n_humans_per_site,
                grid = grid),
           "../data/power_grid.json", auto_unbox = TRUE, pretty = TRUE, digits = 4)
cat("\nWrote power_grid.json\n")
