#' @title Stochastic SIRS simulator (single population)
#'
#' @description
#' Simulate a single-population SIRS process for `n_times` days, tracking S/I/R
#' **counts** internally (for binomial transitions) and returning **proportions**
#' plus daily **cases** (counts). Transmission `beta` may be constant or
#' time-varying (length `n_times`). Optional reporting (`alpha`) thins cases.
#'
#' @param n_times Integer (\eqn{\ge} 2). Number of days.
#' @param pop Integer/numeric (\eqn{>} 0). Total population size (closed).
#' @param I_init Integer/numeric (\eqn{\ge} 0, \eqn{\le} `pop`). Initial infected count.
#' @param beta Numeric scalar or vector (`n_times`). Transmission rate(s) per day
#'   (must be finite and \eqn{\ge} 0).
#' @param gamma Numeric in \[0,1]. Daily recovery probability.
#' @param omega Numeric in \[0,1]. Daily waning probability (R→S).
#' @param epsilon Numeric (\eqn{\ge} 0). External infection pressure.
#' @param alpha `NULL` or numeric in \[0,1]. Reporting probability for thinning cases.
#' @param n_sims Integer (\eqn{\ge} 1). Number of independent simulation runs (columns).
#' @param seed Optional integer. RNG seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{time}{Integer vector `1:n_times`.}
#'   \item{proportions}{Numeric array `[time x sims x state]` with `state` ∈ {`"S"`, `"I"`, `"R"`}.}
#'   \item{cases}{Integer matrix `[time x sims]` of daily (reported) cases.}
#'   \item{params}{List echoing inputs (including the full `beta` vector).}
#' }
#'
#' @examples
#' # 50 runs with constant beta
#' out <- simulate_sirs_stoch(n_times = 200, pop = 1e5, I_init = 20,
#'                      beta = 0.16, gamma = 1/7, omega = 1/30,
#'                      epsilon = 1e-4, n_sims = 50, seed = 42)
#' dim(out$proportions)  # 200 x 50 x 3
#'
#' # Time-varying beta (e.g., seasonal) and reported cases at 60%
#' # b <- make_beta(365, mode = "seasonal", base = 0.18, amplitude = 0.25, phase = 30)
#' # out2 <- simulate_sirs_stoch(365, 2e5, 15, beta = b, gamma = 1/7, omega = 1/60,
#' #                       epsilon = 1e-4, alpha = 0.6, n_sims = 25, seed = 1)
#'
#' @export
simulate_sirs_stoch <- function(
    n_times = 365,     # total number of time steps (e.g., days)
    pop     = 100000,  # population size (closed system)
    I_init  = 10,      # initial infected COUNT at day 1
    beta    = 0.16,    # scalar (constant) OR vector length n_times (time-varying)
    gamma   = 1/7,     # recovery prob per day (e.g., 1/7 -> avg 7 days infectious)
    omega   = 1/30,    # waning prob per day from R->S (e.g., 1/30 -> avg 30 days immune)
    epsilon = 0,       # external infection pressure added to the hazard
    alpha   = NULL,    # optional reporting probability (thins cases)
    n_sims  = 1,       # number of independent simulation runs (columns)
    seed    = NULL     # optional RNG seed for reproducibility
) {
  # ---- checks ----
  stopifnot(n_times >= 2, pop > 0, I_init >= 0, I_init <= pop, n_sims >= 1)
  stopifnot(is.numeric(gamma) && gamma >= 0 && gamma <= 1)
  stopifnot(is.numeric(omega) && omega >= 0 && omega <= 1)
  stopifnot(is.numeric(epsilon) && epsilon >= 0)
  if (!is.null(alpha)) stopifnot(is.numeric(alpha) && alpha >= 0 && alpha <= 1)
  
  # --- make beta a full-length vector ---
  if (length(beta) == 1) {
    beta_vec <- rep(beta, n_times)
  } else {
    if (length(beta) != n_times) stop("beta must be length 1 or n_times.")
    beta_vec <- beta
  }
  if (any(!is.finite(beta_vec)) || any(beta_vec < 0))
    stop("beta must be finite and non-negative.")
  
  # Set random seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # ---- state arrays as COUNTS (dimensions: time x sims) ----
  S <- matrix(0L, n_times, n_sims)
  I <- matrix(0L, n_times, n_sims)
  R <- matrix(0L, n_times, n_sims)
  cases <- matrix(0L, n_times, n_sims)
  
  # Initial conditions at t = 1
  S[1,] <- pop - I_init
  I[1,] <- I_init
  R[1,] <- 0L
  
  # Helper: clamp probs to [0,1]
  pr <- function(x) pmax(pmin(x, 1), 0)
  
  # ---- main loop (always stochastic) ----
  for (t in 2:n_times) {
    # Per-day infection probability from hazard:
    # lambda = 1 - exp(-(beta * I/pop + epsilon))
    lambda  <- 1 - exp(-(beta_vec[t-1] * I[t-1,] / pop + epsilon))
    
    # Binomial transitions
    new_inf <- rbinom(n_sims, size = S[t-1,], prob = pr(lambda))
    new_rec <- rbinom(n_sims, size = I[t-1,], prob = pr(gamma))
    loss_im <- rbinom(n_sims, size = R[t-1,], prob = pr(omega))
    
    # Update compartments
    S[t,] <- S[t-1,] - new_inf + loss_im
    I[t,] <- I[t-1,] + new_inf - new_rec
    R[t,] <- R[t-1,] - loss_im + new_rec
    
    # Record cases (counts), optionally thinned by alpha
    if (is.null(alpha)) {
      cases[t,] <- new_inf
    } else {
      cases[t,] <- rbinom(n_sims, size = new_inf, prob = pr(alpha))
    }
    
    # Guard rails
    S[t,] <- pmax(S[t,], 0L)
    I[t,] <- pmax(I[t,], 0L)
    R[t,] <- pmax(R[t,], 0L)
  }
  
  # ---- return proportions for S/I/R ----
  props <- array(NA_real_, dim = c(n_times, n_sims, 3),
                 dimnames = list(NULL, NULL, c("S","I","R")))
  props[, , "S"] <- S / pop
  props[, , "I"] <- I / pop
  props[, , "R"] <- R / pop
  # ---- return everything tidy ----
  
  list(
    time        = seq_len(n_times),
    proportions = props,
    cases       = cases,
    params      = list(
      n_times = n_times, pop = pop, I_init = I_init, beta = beta_vec,
      gamma = gamma, omega = omega, epsilon = epsilon, alpha = alpha,
      n_sims = n_sims, seed = seed
    )
  )
}
