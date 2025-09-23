#' adjust_beta
#' @title Scale beta within specified day windows (simple intervention/scenario)
#' @description
#' Multiplies beta by a constant `multiplier` in given day ranges.
#' Accepts a vector `[time]` or matrix `[time x groups]` and preserves the shape.
#' @param beta Numeric vector or matrix of beta values.
#' @param windows data.frame with integer columns `start`, `end` (1-based, inclusive).
#' @param multiplier Numeric scalar (e.g., 0.7 means 30% reduction).
#' @return A beta vector/matrix with the same shape as the input.
#' @examples
#' # B <- adjust_beta(rep(0.2, 365), data.frame(start=50, end=80), 0.7)

# Project: Kids Research Institute — SIRS modelling
# Script: R/adjust_beta.R
# Purpose: Scale β(t) within specified day windows for scenarios/interventions
# Inputs: beta vector or [time x groups] matrix; windows data.frame(start,end)
# Outputs: Beta with same shape, scaled in given windows

adjust_beta <- function(beta, windows, multiplier = 0.8) {
  if (is.vector(beta)) {
    B <- matrix(beta, ncol = 1)
    vec <- TRUE
  } else if (is.matrix(beta)) {
    B <- beta
    vec <- FALSE
  } else stop("beta must be a vector or matrix.")
  
  if (!is.data.frame(windows) || !all(c("start", "end") %in% names(windows))) {
    stop("windows must be a data.frame with columns: start, end")
  }
  
  for (i in seq_len(nrow(windows))) {
    s <- max(1, as.integer(windows$start[i]))
    e <- min(nrow(B), as.integer(windows$end[i]))
    if (s <= e) B[s:e, ] <- B[s:e, ] * multiplier
  }
  
  if (vec) drop(B[, 1]) else B
}
