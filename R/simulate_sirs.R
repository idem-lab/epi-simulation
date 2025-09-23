#' simulate_sirs: Discrete-time SIRS simulator (proportions in S, I, R)
#'
#' @param n_times Integer, number of days to simulate (>= 2).
#' @param pop Integer, total closed population (> 0).
#' @param I_init Integer, initial infected count (>= 0, <= pop).
#' @param omega Numeric, waning rate R->S per day (>= 0).
#' @param gamma Numeric, recovery rate I->R per day (>= 0).
#' @param beta_vec Numeric vector length n_times, contact*transmission per day.
#' @param seed Optional integer for reproducibility.
#' @param shock_fun Optional function(t, S_prev, I_prev, R_prev) -> list(S=..., I=..., R=...).
#'        If provided, it is called at each t BEFORE the update; return NULL to do nothing.
#'
#' @return data.frame with columns: t, S, I, R (proportions), N_t (incident count), beta.

# Project: Kids Research Institute — SIRS modelling
# Script: R/simulate_sirs.R
# Purpose: Discrete-time SIRS simulator (single population)
# Inputs: n_times, pop, I_init, omega, gamma, beta_vec, seed, shock_fun
# Outputs: data.frame: t, S, I, R (props), N_t (counts), beta

simulate_sirs <- function(
    n_times,
    pop,
    I_init,
    omega,
    gamma,
    beta_vec,
    seed = NULL,
    shock_fun = NULL
) {
  stopifnot(n_times >= 2, pop > 0, I_init >= 0, I_init <= pop, omega >= 0, gamma >= 0)
  if (is.null(beta_vec)) stop("beta_vec is required and cannot be NULL.")
  stopifnot(length(beta_vec) == n_times, is.numeric(beta_vec), all(beta_vec >= 0))
  
  if (!is.null(seed)) set.seed(seed)
  
  # compartments stored as proportions; incidence as counts
  S <- numeric(n_times)
  I <- numeric(n_times)
  R <- numeric(n_times)
  N_t <- numeric(n_times)
  
  S[1]  <- 1 - I_init / pop
  I[1]  <- I_init / pop
  R[1]  <- 0
  N_t[1] <- I[1] * pop
  
  for (t in 2:n_times) {
    # optional exogenous shock applied to state at t-1 (before advancing to t)
    if (!is.null(shock_fun)) {
      adj <- shock_fun(t, S[t - 1], I[t - 1], R[t - 1])
      if (is.list(adj) && all(c("S", "I", "R") %in% names(adj))) {
        S[t - 1] <- adj$S
        I[t - 1] <- adj$I
        R[t - 1] <- adj$R
      }
    }
    
    beta_prev <- beta_vec[t - 1]
    
    # discrete SIRS updates
    S[t] <- S[t - 1] - beta_prev * S[t - 1] * I[t - 1] + omega * R[t - 1]
    I[t] <- I[t - 1] + beta_prev * S[t - 1] * I[t - 1] - gamma * I[t - 1]
    R[t] <- R[t - 1] - omega * R[t - 1] + gamma * I[t - 1]
    
    # incident infections (emergent)
    N_t[t] <- pop * beta_prev * S[t - 1] * I[t - 1]
  }
  
  data.frame(
    t    = seq_len(n_times),
    S = S, I = I, R = R,
    N_t  = N_t,
    beta = beta_vec
  )
}
