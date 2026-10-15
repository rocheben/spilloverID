
dir.create("../../../data/suppl_S5_uncertainty", showWarnings = FALSE, recursive = TRUE)
# Probability that each host's Pi exceeds the published (0.55) and the
# recalibrated application-regime threshold, under the corrected bootstrap.
boot <- readRDS("../../../data/suppl_S5_uncertainty/boot_raw.rds")
thr  <- read.csv("../../../data/suppl_S5_uncertainty/threshold_recalibration.csv")
t_pub  <- 0.55
t_rec  <- thr$equivalent_threshold[thr$regime == "reconstruction11" &
                                   thr$published_threshold == 0.55]
CASE_TRUE <- c(sn = "P. maniculatus", nv = "Pteropus", rv = "Wild ungulates")
rows <- list()
for (key in names(boot$boots)) {
  p <- strsplit(key, "|", fixed = TRUE)[[1]]; nm <- p[1]; mo <- p[2]
  bm <- boot$boots[[key]]; bm <- bm[rowSums(is.na(bm)) == 0, , drop = FALSE]
  for (k in seq_len(ncol(bm))) rows[[length(rows)+1]] <- data.frame(
    case = nm, host = colnames(bm)[k], scheme = mo,
    is_true_host = colnames(bm)[k] == unname(CASE_TRUE[nm]),
    pi_point = as.numeric(boot$point[[nm]]$pi[k]),
    p_exceed_055 = mean(bm[, k] >= t_pub),
    p_exceed_recalibrated = mean(bm[, k] >= t_rec),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
out$threshold_recalibrated <- round(t_rec, 4)
write.csv(out, "../../../data/suppl_S5_uncertainty/case_threshold_exceedance.csv", row.names = FALSE)
cat("recalibrated threshold =", round(t_rec, 4), "\n\n")
print(out[out$scheme %in% c("published","corrected"), ], row.names = FALSE, digits = 3)
