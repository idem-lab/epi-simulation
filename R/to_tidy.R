#' to_tidy
#' @title Convert simulation output to a tidy long table
#' @description
#' Harmonizes different sim outputs into a single long format with columns:
#' `time, group, sim, state (S/I/R/incidence), value`.
#' @param sim A simulation output (deterministic single-pop, stochastic multi-run, or multi-pop).
#' @return data.frame in long (tidy) format.
#' @examples
#' # df <- to_tidy(sim); head(df)

# Project: Kids Research Institute — SIRS modelling
# Script: R/to_tidy.R
# Purpose: Convert sim outputs to one tidy long table
# Inputs: sim (deterministic | stochastic | multi-pop)
# Outputs: data.frame: time, group, sim, state(S/I/R/incidence), value

to_tidy <- function(sim) {
  pack <- function(t, g, s, S, I, R, inc) {
    rbind(
      data.frame(time = t, group = g, sim = s, state = "S", value = S),
      data.frame(time = t, group = g, sim = s, state = "I", value = I),
      data.frame(time = t, group = g, sim = s, state = "R", value = R),
      data.frame(time = t, group = g, sim = s, state = "incidence", value = inc)
    )
  }
  
  # Deterministic single-pop (vectors S/I/R)
  if (!is.null(sim$S) && is.vector(sim$S) && is.null(sim$proportions)) {
    return(pack(sim$time, 1, 1, sim$S, sim$I, sim$R, sim$incidence))
  }
  
  # Stochastic multi-run (array proportions[time, sim, ])
  if (!is.null(sim$proportions)) {
    M <- dim(sim$proportions)[2]
    out <- lapply(seq_len(M), function(j) {
      pack(sim$time, 1, j,
           S   = sim$proportions[, j, "S"],
           I   = sim$proportions[, j, "I"],
           R   = sim$proportions[, j, "R"],
           inc = sim$cases[, j])
    })
    return(do.call(rbind, out))
  }
  
  # Multi-pop deterministic (matrices S/I/R/incidence)
  if (is.matrix(sim$S) && is.matrix(sim$inci