#' ============================================================================
#' Bootstrap convergence check (Supplementary Figure 7)
#'
#' Width of the 95% site-level bootstrap percentile interval on the true host's
#' composite score as the number of bootstrap replicates B increases, for the
#' three reconstructed outbreaks. A single run of B_MAX replicates is drawn per
#' case and the interval is recomputed on nested prefixes of that run, so the
#' curve shows the convergence of one bootstrap sequence rather than the
#' between-run variability of independent sequences.
#'
#' Requires ../data/case_results.rds (produced by run_cases.R).
#' Writes ../data/boot_converge.json
#' Runtime: ~20-30 min.
#' ============================================================================

suppressPackageStartupMessages({
  library(spilloverID)
  library(jsonlite)
})

set.seed(20260524)

B_MAX  <- 500
B_GRID <- c(25, 50, 100, 150, 200, 300, 400, 500)

cases <- readRDS("../data/case_results.rds")

out <- list()
for (key in c("sn", "nv", "rv")) {
  cc    <- cases[[key]]
  hosts <- names(cc$result$pi)
  ti    <- which(hosts == cc$true_host)
  cat("----", cc$case, "----\n")
  res <- suppressWarnings(compute_pi(cc$data, bootstrap = B_MAX))
  bm  <- res$bootstrap[, ti]
  widths <- sapply(B_GRID, function(B) {
    q <- stats::quantile(bm[seq_len(B)], c(0.025, 0.975), na.rm = TRUE)
    as.numeric(q[2] - q[1])
  })
  for (i in seq_along(B_GRID))
    cat(sprintf("  B=%3d  width=%.4f\n", B_GRID[i], widths[i]))
  out[[key]] <- list(case = cc$case, true_host = cc$true_host,
                     B = B_GRID, width = widths,
                     width_at_Bmax = widths[length(widths)])
}

write_json(out, "../data/boot_converge.json",
           auto_unbox = TRUE, pretty = TRUE, digits = 4)
cat("\nWrote boot_converge.json\n")
