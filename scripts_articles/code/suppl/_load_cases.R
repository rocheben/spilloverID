#' Shared loader for the supplementary analyses.
#'
#' Loads the installed spilloverID package and the three outbreak-reconstruction
#' constructors defined in code/run_cases.R, extracted verbatim so that the
#' supplementary analyses build exactly the datasets the main analysis scores.
suppressPackageStartupMessages(library(spilloverID))

.find_code <- function() {
  for (d in c("..", "../..", "../../..", "."))
    if (file.exists(file.path(d, "run_cases.R"))) return(d)
  stop("run_cases.R not found: run this script from its own directory")
}
.CODE_DIR <- .find_code()

.src  <- readLines(file.path(.CODE_DIR, "run_cases.R"))
.from <- grep("^# ---- Dataset constructor", .src)[1]
.to   <- grep("^# ---- Run all three cases", .src)[1] - 1L
.defs <- tempfile(fileext = ".R")
writeLines(.src[.from:.to], .defs)
source(.defs)

build_cases <- function(seed = 20260524) {
  set.seed(seed)                       # same seed and same order as run_cases.R
  list(sn = case_sin_nombre(), nv = case_nipah(), rv = case_rvf())
}

CASE_TRUE <- c(sn = "P. maniculatus", nv = "Pteropus", rv = "Wild ungulates")
CASE_LAB  <- c(sn = "Sin Nombre (USA, 1993)",
               nv = "Nipah (Malaysia, 1998-99)",
               rv = "Rift Valley fever (Kenya, 2006-07)")

CRITF <- list(
  C1 = criterion_C1_spatial, C2 = criterion_C2_reservoir_signatures,
  C3 = criterion_C3_individual, C4 = criterion_C4_dose_response,
  C5 = criterion_C5_cryptic_fraction, C6 = criterion_C6_interaction,
  C7 = criterion_C7_foi_coherence, C8 = criterion_C8_temporal_precedence,
  C9 = criterion_C9_mediation, C10 = criterion_C10_mechanistic,
  C11 = criterion_C11_specificity, C12 = criterion_C12_inter_epidemic,
  C13 = criterion_C13_endemic_age, C14 = criterion_C14_spatial_independence
)
