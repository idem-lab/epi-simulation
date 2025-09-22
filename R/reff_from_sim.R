#' reff_from_sim
#' @title Compute effective reproduction number R_eff(t)
#' @description
#' Works with three shapes of simulation outputs produced in this repo:
#' - Deterministic single-pop (`simulate_sirs_det`): list with vectors `S`, `I`, `R`,
#'   `incidence`, and `params$beta`, `params$gamma`.
#' - Stochastic multi-run (`simulate_sirs`): list with array `proportions[time, sim, state]`,
#'   matrix `cases[time, sim]`, and `params$beta`, `params$gamma`.
#' - Multi-pop deterministic (`simulate_sirs_multi`): list with matrices `S`, `I`, `R`,
#'   `incidence`, and `params$beta` (vector or [time x groups]) and `params$gamma`.
#'
#' Uses the simple approximation R_eff(t) ≈ S(t) * beta(t) / gamma.
#' For multi-pop outputs this is computed per group and **ignores** cross-mixing (C),
#' which is acceptable for quick diagnostics.
#'
#' @param sim A simulation output list as described above.
#' @return A data.frame with columns:
#'   - deterministic single-pop: `time, Reff`
#'   - stochastic:                `time, sim, Reff`
#'   - multi-pop:                 `time, group, Reff`
#' @examples
#' # df <- reff_from_sim(sim); head(df)

# Project: Kids Research Institute — SIRS modelling
# Script: R/reff_from_sim.R
# Purpose: Compute R_eff(t) for deterministic, stochastic, and multi-pop outputs
# Inputs: sim list with params$beta and params$gamma; uses S (or proportions[,, "S"])
# Outputs: data.frame(time, Reff) or with sim/group columns

reff_from_sim <- function(sim) {
  # helpers to pull parameters
  get_param <- function(x, k) {
    if (is.list(x$params) && !is.null(x$params[[k]])) return(x$params[[k]])
    stop(k, " not found in sim$params")
  }
  gamma <- get_param(sim, "gamma")
  beta  <- get_param(sim, "beta")
  
  # Case A: deterministic single-pop (vectors S/I/R)
  if (!is.null(sim$S) && is.vector(sim$S) && is.null(sim$proportions)) {
    S <- sim$S
    if (length(beta) == 1) beta <- rep(beta, length(S))
    if (length(beta) != length(S)) stop("beta length must match the simulation length.")
    Reff <- S * (beta / gamma)
    return(data.frame(time = sim$time, Reff = Reff))
  }
  
  # Case B: stochastic multi-run (array proportions[time, sim, ])
  if (!is.null(sim$proportions)) {
    T <- dim(sim$proportions)[1]
    M <- dim(sim$proportions)[2]
    if (length(beta) == 1) beta <- rep(beta, T)
    if (length(beta) != T) stop("beta length must match time for stochastic sim.")
    out <- lapply(seq_len(M), function(j) {
      S <- sim$proportions[, j, "S"]
      data.frame(time = sim$time, sim = j, Reff = S * (beta / gamma))
    })
    return(do.call(rbind, out))
  }
  
  # Case C: multi-pop deterministic (matrices S/I/R)
  if (is.matrix(sim$S)) {
    S <- sim$S
    T <- nrow(S); G <- ncol(S)
    if (is.vector(beta)) beta <- matrix(beta, nrow = T, ncol = G)
    if (!is.matrix(beta) || nrow(beta) != T || ncol(