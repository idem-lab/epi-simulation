# =============================================================================
# Project: Kids Research Institute — SIRS modelling
# File   : R/cum_metrics.R
# Purpose: Cumulative incidence and attack-rate helpers
# Inputs : sim (deterministic | stochastic | multi-pop; single or multi-run)
# Outputs:
#   - cumulative_incidence(sim) -> data.frame with time and cum_* plus
#     optional columns: sim (simulation index) and/or group (group index)
#   - attack_rate(sim)          -> scalar / vector / matrix (see below)
# Notes : Pure base R (no dplyr). ASCII only.
# =============================================================================

# Small helper: x %||% y  -> return x unless it's NULL
`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal: extract population vector per group (or scalar for single-pop)
.get_pop_vec <- function(sim) {
  # Prefer explicit vector if present
  pv <- sim$params$pop_vec %||% sim$pop_vec
  if (!is.null(pv)) return(as.numeric(pv))
  
  # Fall back to scalar pop in params or top-level
  p  <- sim$params$pop %||% sim$pop
  if (!is.null(p)) return(as.numeric(p))
  
  # Last resort: 1 (prevents division by zero)
  1
}

# -----------------------------------------------------------------------------
# cumulative_incidence()
# -----------------------------------------------------------------------------
# Returns cumulative sums of daily infections (preferred: 'incidence').
# If 'incidence' is absent, falls back to 'cases' (reported).
# Shapes handled:
# - deterministic single-pop: vector [time]
# - deterministic multi-pop : matrix [time x groups]
# - stochastic single-pop   : matrix [time x sims]
# - stochastic multi-pop    : array  [time x sims x groups]
#
# Output is always a data.frame with columns:
# - time
# - cum_incidence  (or cum_cases if incidence absent)
# - optional: sim (for stochastic), group (for multi-pop)
# -----------------------------------------------------------------------------
cumulative_incidence <- function(sim) {
  use_field <- if (!is.null(sim$incidence)) "incidence" else if (!is.null(sim$cases)) "cases" else NULL
  if (is.null(use_field)) stop("No 'incidence' or 'cases' found in sim for cumulative_incidence().")
  
  x <- sim[[use_field]]
  nm <- if (use_field == "incidence") "cum_incidence" else "cum_cases"
  
  # 1) Vector [time]
  if (is.vector(x) && is.null(dim(x))) {
    return(data.frame(
      time = sim$time,
      ..tmp.. = cumsum(as.numeric(x))
      , check.names = FALSE, fix.empty.names = FALSE,
      row.names = NULL, stringsAsFactors = FALSE)[, c("time", "..tmp.."),
                                