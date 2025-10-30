#' @title Deterministic multi-population SIRS simulator
#' @description
#' Simulate a deterministic SIRS model across \eqn{P} subpopulations over
#' `n_times` discrete steps. Returns S, I, R **proportions** and daily
#' **incidence (counts)** per group. Transmission can vary over time and/or
#' by group. **No contact matrix is used**; each group’s force of infection is
#' within-group only: \eqn{\lambda_p(t) = \beta_{p}(t)\, I_p(t)}.
#'
#' @param n_times Integer (\eqn{\ge} 2). Number of time steps (days).
#' @param pop_vec Numeric vector (length `P`, all \eqn{>} 0). Population sizes per group.
#' @param I_init Numeric/integer vector (length `P`, each in \[0, `pop_vec[p]`]).
#'   Initial infected **counts** per group.
#' @param beta_mat Numeric scalar, vector (`n_times`), or matrix (`[n_times x P]`)
#'   giving transmission rates.
#' @param gamma Numeric in \[0,1]. Recovery rate per day (common across groups).
#' @param omega Numeric in \[0,1]. Waning rate per day from R to S (common).
#' @param seed Optional integer. Only used to fix upstream stochastic generation;
#'   this simulator itself is deterministic.
#'
#' @examples
#' # Two groups, constant beta, no mixing between groups (decoupled)
#' out <- simulate_sirs_multi(
#'   n_times = 90,
#'   pop_vec = c(5e4, 5e4),
#'   I_init  = c(10, 10),
#'   beta_mat = 0.16,
#'   gamma = 1/7, omega = 1/30
#' )
#' str(out$S); str(out$incidence)
#'
#' # Time-varying, group-specific beta
#' # b1 <- make_beta(120, mode = "seasonal", base = 0.18, amplitude = 0.25, phase = 30)
#' # b2 <- make_beta(120, mode = "constant", value = 0.16)
#' # B  <- cbind(b1, b2)  # [n_times x P]
#' # out2 <- simulate_sirs_multi(120, c(8e4, 6e4), c(8, 6), B, 1/7, 1/60)
#'
#' @export
simulate_sirs_multi <- function(
    n_times   = 365,                 # number of time steps (e.g., days)
    pop_vec   = c(50000, 50000),     # population sizes per group (length = P)
    I_init    = c(10, 10),           # initial infected COUNTS per group (length = P)
    beta_mat  = NULL,                # transmissibility: scalar / vector(n_times) / matrix[n_times x P]
    gamma     = 1/7,                 # recovery rate per day (same for all groups here)
    omega     = 1/30,                # waning rate per day (R -> S, same for all groups)
    seed      = 123                  # RNG seed for reproducibility (no randomness here, but ok)
) {
  # ----- basic sizes -----
  P <- length(pop_vec)                              # number of groups/subpopulations
  stopifnot(length(I_init) == P)                    # ensure I_init matches pop_vec
  stopifnot(all(I_init >= 0 & I_init <= pop_vec))   # valid initial infected counts
  stopifnot(n_times >= 2, P >= 1)
  
  if (!is.null(seed)) set.seed(seed)
  
  # ----- handle beta_mat options -----
  # Accept:
  #  - scalar: same value for all days & groups
  #  - vector length n_times: same across groups but time-varying
  #  - matrix [n_times x P]: time-varying AND group-specific
  if (is.null(beta_mat)) {
    beta_mat <- matrix(0.16, nrow = n_times, ncol = P)  # default constant 0.16
  } else if (length(beta_mat) == 1) {
    beta_mat <- matrix(beta_mat, nrow = n_times, ncol = P)
  } else if (is.vector(beta_mat) && length(beta_mat) == n_times) {
    beta_mat <- matrix(beta_mat, nrow = n_times, ncol = P)
  } else if (is.matrix(beta_mat)) {
    stopifnot(nrow(beta_mat) == n_times, ncol(beta_mat) == P)
  } else {
    stop("beta_mat must be scalar, vector of length n_times, or matrix [n_times x P].")
  }
  if (any(!is.finite(beta_mat)) || any(beta_mat < 0))
    stop("beta_mat must be finite and non-negative.")
  
  # ----- initialize states as PROPORTIONS (one column per group) -----
  S <- matrix(0, nrow = n_times, ncol = P)   # susceptible proportion in each group
  I <- matrix(0, nrow = n_times, ncol = P)   # infected proportion in each group
  R <- matrix(0, nrow = n_times, ncol = P)   # recovered proportion in each group
  
  for (p in 1:P) {
    S[1, p] <- 1 - I_init[p] / pop_vec[p]    # convert initial infected counts to proportions
    I[1, p] <- I_init[p] / pop_vec[p]
    R[1, p] <- 0
  }
  
  # Daily incidence (COUNTS) per group
  incidence <- matrix(0, nrow = n_times, ncol = P)
  incidence[1, ] <- I[1, ] * pop_vec         # day-1 = seeded infections (optional convention)
  
  # ----- main loop over time -----
  for (t in 2:n_times) {
    # Force of infection λ_p(t-1) = β_{p}(t-1) * I_p(t-1)
    lambda <- beta_mat[t-1, ] * I[t-1, ]  # length P, element-wise
    
    for (p in 1:P) {
      # Flows for group p at day t-1 -> t
      inf_flow   <- lambda[p] * S[t-1, p]   # new infections (proportion)
      recov_flow <- gamma    * I[t-1, p]    # recoveries (proportion)
      wane_flow  <- omega    * R[t-1, p]    # loss of immunity (proportion)
      
      # Euler updates (discrete-time SIRS)
      S[t, p] <- S[t-1, p] - inf_flow + wane_flow
      I[t, p] <- I[t-1, p] + inf_flow - recov_flow
      R[t, p] <- R[t-1, p] + recov_flow - wane_flow
      
      # Numerical hygiene: clamp and renormalize to keep S+I+R ≈ 1
      S[t, p] <- max(S[t, p], 0)
      I[t, p] <- max(I[t, p], 0)
      R[t, p] <- max(R[t, p], 0)
      tot <- S[t, p] + I[t, p] + R[t, p]
      if (!isTRUE(all.equal(tot, 1, tolerance = 1e-10))) {
        S[t, p] <- S[t, p] / tot
        I[t, p] <- I[t, p] / tot
        R[t, p] <- R[t, p] / tot
      }
      
      # Incidence as COUNTS for group p on day t
      incidence[t, p] <- pop_vec[p] * inf_flow
    }
  }
  
  # ----- return a tidy list -----
  list(
    time = seq_len(n_times),     # 1..n_times
    S = S, I = I, R = R,         # proportions per group over time
    incidence = incidence,       # counts per day per group
    params = list(
      n_times = n_times,
      pop_vec = pop_vec,
      I_init  = I_init,
      gamma   = gamma,
      omega   = omega,
      beta    = beta_mat
    )
  )
}
