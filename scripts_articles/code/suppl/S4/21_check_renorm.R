# Does the discount renormalize the weight vector or not? Both conventions are
# defensible; they give different q*. Report which one reproduces the audit.
RI <- readRDS("reliability_inputs.rds")
source("../_load_cases.R")
CASES <- lapply(names(CASE_LAB), function(k) list(label = unname(CASE_LAB[k]), data = NULL, true_host = unname(CASE_TRUE[k])))
names(CASES) <- names(CASE_LAB)
.blt <- build_cases()
for (k in names(CASES)) CASES[[k]]$data <- .blt[[k]]
BASE <- lapply(CASES, function(cs) suppressWarnings(compute_pi(cs$data, bootstrap = 0)))
CRIT <- lapply(BASE, function(b) b$criteria)
DEP <- CRIT
W0 <- default_weights(); CRITNM <- names(W0)
TRUTH <- c(sn = "P. maniculatus", nv = "Pteropus", rv = "Wild ungulates")
QG <- seq(0, 1, by = 0.0005)

margin_at <- function(cm, w, renorm) {
  keep <- !apply(cm, 2, function(x) all(is.na(x)))
  w[!keep] <- 0
  if (renorm) w <- w / sum(w)
  cs <- cm; cs[is.na(cs)] <- 0
  pv <- stats::setNames(as.numeric(cs %*% w[CRITNM]), rownames(cm))
  pv
}
tipq <- function(cm, truth, target, renorm) {
  m <- sapply(QG, function(q) {
    w <- W0; w[target] <- w[target] * q
    pv <- margin_at(cm, w, renorm)
    comp <- names(sort(pv, decreasing = TRUE))[names(sort(pv, decreasing = TRUE)) != truth][1]
    pv[truth] - pv[comp]
  })
  i <- which(rev(m) <= 0)[1]
  if (is.na(i)) NA_real_ else rev(QG)[i]
}
res <- expand.grid(case = c("sn","nv","rv"), target = c("triad","C12"),
                   renorm = c(TRUE, FALSE), stringsAsFactors = FALSE)
res$q_star <- mapply(function(cs, tg, rn)
  tipq(DEP[[cs]], TRUTH[[cs]], if (tg == "triad") c("C12","C13","C14") else "C12", rn),
  res$case, res$target, res$renorm)
print(res, row.names = FALSE)
