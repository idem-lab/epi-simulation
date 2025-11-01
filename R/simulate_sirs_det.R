#' @title Deterministic SIRS simulator (single population)
#'
#' @description
#' Simulate a deterministic SIRS model (no births/deaths/migration) over
#' `n_times` discrete time steps, returning S, I, R **proportions** and
#' daily **incidence (integer counts)**. The transmission rate `beta` can be a
#' scalar (constant) or a vector of length `n_times` (time-varying).
#'
#' @details
#' This is a deterministic, continuous-flow SIRS model. Internally, new infections
#' on day *t* are computed as an expected value \eqn{\text{pop} \times \beta_t S_{t-1} I_{t-1}}.
#' For convenience, the returned `incidence` is the integer-rounded version of that
#' expected count (`as.integer(round(...))`) for reporting/plotting. If you need
#' exact expected (non-integer) counts, compute them separately from the returned
#' S, I, R and parameters, or use a stochastic simulator.
#'
#' @param n_times Integer (≥ 2). Number of time steps (days).
#' @param pop Integer/numeric (> 0). Total population size.
#' @param I_init Integer/numeric (≥ 0). Initial infected **count**.
#' @param beta Numeric. Either a scalar (constant) or a vector of length `n_times`
#'   giving the transmission rate each day (must be finite and ≥ 0).
#' @param gamma Numeric in [0,1]. Recovery rate per day (e.g., `1/7`).
#' @param omega Numeric in [0,1]. Waning rate per day from R to S (e.g., `1/30`).
#' @param seed Optional integer. If supplied, used to set the RNG seed.
#'
#' @return A list with components:
#' \describe{
#'   \item{time}{Integer vector `1:n_times`.}
#'   \item{S, I, R}{Numeric vectors of length `n_times` (proportions).}
#'   \item{incidence}{**Integer** vector of length `n_times` (daily counts of new infections,
#'         computed as `round(pop * beta * S * I)`).}
#'   \item{params}{List of inputs used (including the full `beta` vector) for reproducibility.}
#' }
#'
#' @examples
#' # Constant-beta run (incidence is integer-rounded for display/reporting)
#' out1 <- simulate_sirs_det(n_times = 60, pop = 1e5, I_init = 10,
#'                           beta = 0.16, gamma = 1/7, omega = 1/30)
#' head(data.frame(day = out1$time, S = out1$S, I = out1$I,
#'                 R = out1$R, inc = out1$incidence))
#'
#' @export

simulate_sirs_det <- function(
    n_times = 365,        # total number of time steps (e.g., days)
    pop     = 100000,     # closed population size
    I_init  = 10,         # initial infected (count)
    beta    = 0.16,       # transmission: scalar OR vector length n_times
    gamma   = 1/7,        # recovery rate per day
    omega   = 1/30,       # waning rate per day (R -> S)
    seed    = 123         # optional RNG seed (useful if beta was generated randomly)
) {
  # --- Checks ---
  stopifnot(n_times >= 2, pop > 0, I_init >= 0, I_init <= pop)
  stopifnot(is.numeric(gamma), gamma >= 0, gamma <= 1)
  stopifnot(is.numeric(omega), omega >= 0, omega <= 1)
  
  # --- Beta as full-length vector ---
  if (length(beta) == 1) {
    beta_vec <- rep(beta, n_times)
  } else {
    if (length(beta) != n_times) stop("beta must be length 1 or n_times.")
    beta_vec <- beta
  }
  if (any(!is.finite(beta_vec)) || any(beta_vec < 0)) stop("beta must be finite and non-negative.")
  
  # Optional seed
  if (!is.null(seed)) set.seed(seed)
  
  # --- Storage for proportions ---
  S <- numeric(n_times)
  I <- numeric(n_times)
  R <- numeric(n_times)
  
  # Initial conditions (day 1)
  S[1] <- 1 - I_init / pop
  I[1] <- I_init / pop
  R[1] <- 0
  
  # --- Daily incidence (expected counts) ---
  incidence <- numeric(n_times)
  incidence[1] <- I[1] * pop  # seeded infections shown on day 1
  
  # --- Iterate days ---
  for (t in 2:n_times) {
    # Flows (as proportions of population)
    inf_flow   <- beta_vec[t-1] * S[t-1] * I[t-1]  # new infections
    recov_flow <- gamma * I[t-1]                   # recoveries
    wane_flow  <- omega * R[t-1]                   # waning immunity
    
    # Update compartments
    S[t] <- S[t-1] - inf_flow + wane_flow
    I[t] <- I[t-1] + inf_flow - recov_flow
    R[t] <- R[t-1] + recov_flow - wane_flow
    
    # Numerical hygiene
    S[t] <- max(S[t], 0); I[t] <- max(I[t], 0); R[t] <- max(R[t], 0)
    tot <- S[t] + I[t] + R[t]
    if (!isTRUE(all.equal(tot, 1, tolerance = 1e-10))) {
      S[t] <- S[t] / tot; I[t] <- I[t] / tot; R[t] <- R[t] / tot
    }
    
    # Expected daily incidence (counts)
    incidence[t] <- pop * inf_flow
  }
  
  # --- Integer copy for presentation ---
  incidence_int <- as.integer(round(incidence, 0))
  
  # --- Return ---
  list(
    time = seq_len(n_times),
    S = S, I = I, R = R,
    incidence = incidence_int,           # continuous expected counts
    params = list(
      n_times = n_times,
      pop     = pop,
      I_init  = I_init,
      gamma   = gamma,
      omega   = omega,
      beta    = beta_vec
    )
  )
}
