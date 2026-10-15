#' ============================================================================
#' Robustness of the interpretation threshold (Supplementary Figure 5)
#'
#' Re-runs the no-spillover negative control over a grid crossing the cryptic
#' background force of infection with the serological specificity, and records
#' in each cell the distribution of the maximum composite score and the
#' false-positive rate at the 0.55 threshold.
#'
#' Writes ../data/fig_threshold_grid.csv and ../data/suppl_threshold_sweep.json
#' ============================================================================
suppressPackageStartupMessages({ library(spilloverID); library(parallel) })
N <- as.integer(Sys.getenv("SWEEP_N_REP", "20"))
.dc <- suppressWarnings(detectCores()); if (!is.finite(.dc)) .dc <- 1L
N_CORE <- max(1L, as.integer(Sys.getenv("SWEEP_N_CORE", as.character(.dc))))
alpha_grid <- c(0.0008, 0.0020, 0.0040, 0.0070)
spec_grid  <- c(0.90, 0.94, 0.97, 0.99)
ckpt <- "../data/threshold_sweep_ckpt.rds"
acc <- if (file.exists(ckpt)) readRDS(ckpt) else list()

for (a in alpha_grid) for (sp in spec_grid) {
  key <- sprintf("a%g_s%g", a, sp)
  if (!is.null(acc[[key]])) { cat("skip", key, "\n"); next }
  got <- mclapply(1:N, function(rep) {
    p <- scenario_params("S11_no_spillover")
    p$alpha_cryptic_bg <- a
    p$serology_specificity <- sp
    seed <- 20260602 + rep*7 + round(a*1e4) + round(sp*100)
    d <- tryCatch(simulate_spillover(p, seed = seed), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    r <- tryCatch(suppressWarnings(compute_pi(d, bootstrap = 0)), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    c(max(r$pi, na.rm = TRUE), mean(d$human$sero, na.rm = TRUE))
  }, mc.cores = N_CORE)
  got <- Filter(function(x) is.numeric(x) && length(x) == 2, got)
  maxpi <- vapply(got, `[`, numeric(1), 1); bgsero <- vapply(got, `[`, numeric(1), 2)
  acc[[key]] <- list(alpha = a, spec = sp, n = length(maxpi),
                     bg_sero = mean(bgsero),
                     maxpi_median = as.numeric(median(maxpi)),
                     maxpi_975 = as.numeric(quantile(maxpi, 0.975)),
                     maxpi_max = max(maxpi),
                     fp_050 = mean(maxpi >= 0.50),
                     fp_055 = mean(maxpi >= 0.55),
                     fp_060 = mean(maxpi >= 0.60),
                     maxpi = maxpi)
  saveRDS(acc, ckpt)
  cat(sprintf("done %s: bg_sero=%.1f%% median_maxPi=%.3f fp055=%.1f%%\n",
              key, 100*acc[[key]]$bg_sero, acc[[key]]$maxpi_median,
              100*acc[[key]]$fp_055))
}

grid <- do.call(rbind, lapply(acc, function(z) data.frame(
  alpha = z$alpha, spec = z$spec, n = z$n, bg_sero = round(z$bg_sero, 4),
  maxpi_median = round(z$maxpi_median, 4), maxpi_975 = round(z$maxpi_975, 4),
  maxpi_max = round(z$maxpi_max, 4), fp_050 = round(z$fp_050, 4),
  fp_055 = round(z$fp_055, 4), fp_060 = round(z$fp_060, 4))))
write.csv(grid, "../data/fig_threshold_grid.csv", row.names = FALSE)
if (requireNamespace("jsonlite", quietly = TRUE))
  jsonlite::write_json(acc, "../data/suppl_threshold_sweep.json",
                       auto_unbox = TRUE, pretty = TRUE, digits = 4)
cat("THRESHOLD_SWEEP DONE. Wrote fig_threshold_grid.csv\n")
