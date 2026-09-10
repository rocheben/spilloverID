#' Build a 'spillover_data' object from field data
#'
#' Combines human serological survey, animal sampling data, and site-level
#' covariates into a standardized object accepted by all criterion functions
#' and \code{\link{compute_pi}}.
#'
#' @param human Data frame of human individuals. Must contain columns:
#'   \code{site} (site ID), \code{ind} (individual ID), \code{age} (numeric),
#'   \code{sero} (binary, 1 = seropositive), and one column per candidate host
#'   named \code{exp_h1}, \code{exp_h2}, ... (binary exposure, declared).
#'   Recommended additional columns: \code{sex}, \code{prof} (occupational risk),
#'   \code{residence} (years at current site), \code{household},
#'   \code{env_declared} (binary environmental exposure indicator).
#' @param animal Data frame of animal observations aggregated at site x host
#'   level. Must contain: \code{site}, \code{host} (integer host index, starting at 1),
#'   \code{n_sampled}, \code{n_pcr} (number PCR-positive),
#'   \code{prev_pcr} (PCR prevalence), \code{n_sero}, \code{prev_sero}.
#' @param animal_indiv (Optional) Individual-level animal data with columns
#'   \code{site}, \code{host}, \code{age}, \code{sero}. Required for criteria C2
#'   (age slope), C7 (FOI catalytic).
#' @param animal_ts (Optional) Animal time series for criterion C8.
#'   Columns: \code{site}, \code{host}, \code{year}, \code{prev_pcr_true}.
#' @param human_ts (Optional) Human incidence time series for criterion C8.
#'   Columns: \code{site}, \code{year}, \code{incidence}.
#' @param site Data frame at site level. Must contain \code{site} and any
#'   environmental covariate (e.g. \code{anthrop}). Per-host columns
#'   \code{prev_pcr_h1}, ..., \code{prev_pcr_hK} and similarly for sero must
#'   be derivable from \code{animal} (auto-filled if missing).
#' @param n_hosts Integer number of candidate host species K.
#' @param host_names Character vector of host names (length K).
#'
#' @return A 'spillover_data' object.
#' @export
spillover_data <- function(human, animal, animal_indiv = NULL,
                           animal_ts = NULL, human_ts = NULL,
                           site = NULL, n_hosts = NULL, host_names = NULL) {
  # Infer K if needed
  if (is.null(n_hosts)) {
    exp_cols <- grep("^exp_h[0-9]+$", names(human), value = TRUE)
    n_hosts <- length(exp_cols)
    if (n_hosts == 0)
      stop("Cannot infer n_hosts: no columns named 'exp_h1', 'exp_h2', ... found in human data.")
  }
  if (is.null(host_names)) host_names <- paste0("h", seq_len(n_hosts))

  # Validate required columns
  req_human <- c("site", "age", "sero", paste0("exp_h", seq_len(n_hosts)))
  miss <- setdiff(req_human, names(human))
  if (length(miss) > 0) stop("Missing columns in 'human': ", paste(miss, collapse = ", "))

  req_animal <- c("site", "host", "n_sampled", "prev_pcr", "prev_sero")
  miss <- setdiff(req_animal, names(animal))
  if (length(miss) > 0) stop("Missing columns in 'animal': ", paste(miss, collapse = ", "))

  # Build site-level prevalence if not provided
  if (is.null(site)) {
    site <- data.frame(site = sort(unique(human$site)))
  }
  for (k in seq_len(n_hosts)) {
    col_pcr <- paste0("prev_pcr_h", k)
    col_sero <- paste0("prev_sero_h", k)
    if (!col_pcr %in% names(site)) {
      a_k <- animal[animal$host == k, c("site", "prev_pcr")]
      site[[col_pcr]] <- a_k$prev_pcr[match(site$site, a_k$site)]
    }
    if (!col_sero %in% names(site)) {
      a_k <- animal[animal$host == k, c("site", "prev_sero")]
      site[[col_sero]] <- a_k$prev_sero[match(site$site, a_k$site)]
    }
  }

  obj <- list(human = human, animal = animal,
              animal_indiv = animal_indiv, animal_ts = animal_ts,
              human_ts = human_ts, site = site,
              n_hosts = n_hosts, host_names = host_names)
  class(obj) <- c("spillover_data", "list")
  obj
}

#' @export
print.spillover_data <- function(x, ...) {
  cat("spillover_data object\n")
  cat("  ", x$n_hosts, "candidate hosts:", paste(x$host_names, collapse = ", "), "\n")
  cat("  ", length(unique(x$human$site)), "sites\n")
  cat("  ", nrow(x$human), "humans sampled (", sum(x$human$sero), "seropositive)\n")
  cat("  ", nrow(x$animal), "site-host animal aggregations\n")
  cat("  Optional data:\n")
  cat("    animal_indiv:", !is.null(x$animal_indiv), "\n")
  cat("    animal_ts   :", !is.null(x$animal_ts), "\n")
  cat("    human_ts    :", !is.null(x$human_ts), "\n")
  invisible(x)
}
