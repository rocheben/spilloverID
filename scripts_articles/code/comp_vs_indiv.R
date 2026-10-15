#' ============================================================================
#' Composite score versus each criterion used alone (Supplementary Figure 4)
#'
#' Derived from the per-replicate criterion matrix written by run_benchmark.R,
#' so the comparison uses the same replicates as the benchmark itself rather
#' than a second, smaller simulation pass. For each replicate and each
#' criterion, the true reservoir is scored top-1 when that criterion alone ranks
#' it strictly first (ties count as a miss, since a tied criterion does not
#' identify anything) and top-2 when it ranks it first or second.
#'
#' Requires ../data/benchmark_criteria.csv (run_benchmark.R).
#' Writes   ../data/fig_compindiv_sim.csv, ../data/fig_compindiv_perscen.csv
#' Runtime: seconds.
#' ============================================================================

crit <- paste0("C", 1:14)
d <- utils::read.csv("../data/benchmark_criteria.csv", stringsAsFactors = FALSE)
d <- d[d$scenario != "S11_no_spillover" & !is.na(d$is_true), ]
d$key <- paste(d$scenario, d$rep, sep = "|")

score_one <- function(v, is_true) {
  if (all(is.na(v))) return(c(NA_real_, NA_real_))
  v2 <- v; v2[is.na(v2)] <- -Inf
  rk <- rank(-v2, ties.method = "min")
  k <- which(is_true)
  c(as.numeric(rk[k] == 1 && sum(v2 == max(v2)) == 1), as.numeric(rk[k] <= 2))
}

keys <- unique(d$key)
rows <- lapply(keys, function(kk) {
  s <- d[d$key == kk, ]
  out <- data.frame(scenario = s$scenario[1], rep = s$rep[1])
  for (cn in crit) {
    r <- score_one(s[[cn]], s$is_true)
    out[[paste0(cn, "_t1")]] <- r[1]; out[[paste0(cn, "_t2")]] <- r[2]
  }
  out
})
big <- do.call(rbind, rows)

bench <- jsonlite::fromJSON("../data/benchmark_results.json")$per_rep
bench <- bench[bench$scenario != "S11_no_spillover", ]
big <- merge(big, bench[, c("scenario", "rep", "top1", "top2")],
             by = c("scenario", "rep"))
names(big)[names(big) == "top1"] <- "comp_t1"
names(big)[names(big) == "top2"] <- "comp_t2"

sim <- data.frame(
  item = c("composite", crit),
  sim_top1 = c(round(100 * mean(big$comp_t1), 1),
           sapply(crit, function(cn) round(100 * mean(big[[paste0(cn, "_t1")]], na.rm = TRUE), 1))),
  sim_top2 = c(round(100 * mean(big$comp_t2), 1),
           sapply(crit, function(cn) round(100 * mean(big[[paste0(cn, "_t2")]], na.rm = TRUE), 1))),
  n_evaluable = c(nrow(big),
           sapply(crit, function(cn) sum(!is.na(big[[paste0(cn, "_t1")]])))),
  row.names = NULL)

pos <- unique(big$scenario)
perscen <- do.call(rbind, lapply(pos, function(sc) {
  s <- big[big$scenario == sc, ]
  vals <- sapply(crit, function(cn) mean(s[[paste0(cn, "_t1")]], na.rm = TRUE) * 100)
  vals[!is.finite(vals)] <- NA
  data.frame(scenario = sc,
             composite_top1 = round(mean(s$comp_t1) * 100, 1),
             best_single_top1 = round(max(vals, na.rm = TRUE), 1),
             best_single = crit[which.max(replace(vals, is.na(vals), -1))])
}))

write.csv(sim, "../data/fig_compindiv_sim.csv", row.names = FALSE)
write.csv(perscen, "../data/fig_compindiv_perscen.csv", row.names = FALSE)
cat("Composite top-1 on the benchmark:", sim$sim_top1[1], "%\n")
cat("Wrote fig_compindiv_sim.csv and fig_compindiv_perscen.csv\n")
