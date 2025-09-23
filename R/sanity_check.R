#' sanity_check: Verify that S + I + R ≈ 1 at each time
#'
#' @param sim_df data.frame from simulate_sirs() with columns S, I, R.
#' @param tol Numeric tolerance for equality (default 1e-12).
#' @return TRUE if within tolerance, otherwise FALSE.

# Project: Kids Research Institute — SIRS modelling
# Script: R/sanity_check.R
# Purpose: Verify S + I + R ≈ 1 at each time step
# Inputs: sim data.frame with S, I, R; tol
# Outputs: TRUE/FALSE

sanity_check <- function(sim_df, tol = 1e-12) {
  stopifnot(all(c("S","I","R") %in% names(sim_df)))
  sums <- sim_df$S + sim_df$I + sim_df$R
  isTRUE(all.equal(rep(1, length(sums)), sums, tolerance = tol))
}
