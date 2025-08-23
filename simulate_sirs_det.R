simulate_sirs_det <- function(
    n_times = 365,        # total number of time steps (e.g., days) to simulate
    pop     = 100000,     # total population size (closed population, no births/deaths)
    I_init  = 10,         # initial number of infected individuals
    beta    = 0.16,       # transmission parameter: scalar (constant) OR vector of length n_times
    gamma   = 1/7,        # recovery rate: 1/7 → avg infectious period ≈ 7 days
    omega   = 1/30,       # waning rate: 1/30 → avg immunity duration ≈ 30 days
    seed    = 123         # optional random seed, only useful if beta was generated randomly
) {
  # --- Check inputs ---
  stopifnot(n_times >= 2, pop > 0, I_init >= 0, I_init <= pop)   # sensible ranges
  stopifnot(is.numeric(gamma), gamma >= 0, gamma <= 1)           # gamma must be between 0 and 1
  stopifnot(is.numeric(omega), omega >= 0, omega <= 1)           # omega must be between 0 and 1

  # --- Make beta a full-length vector ---
  # If user gave just one number (constant beta), repeat it for all days.
  # If user gave a vector, it must be exactly length = n_times.
  if (length(beta) == 1) {
    beta_vec <- rep(beta, n_times)
  } else {
    if (length(beta) != n_times) stop("beta must be length 1 or n_times.")
    beta_vec <- beta
  }
  # Check beta values are valid
  if (any(!is.finite(beta_vec)) || any(beta_vec < 0)) stop("beta must be finite and non-negative.")

  # Set random seed if provided (useful when beta was generated with randomness)
  if (!is.null(seed)) set.seed(seed)

  # --- Allocate storage for S, I, R over time ---
  # These are proportions of the total population (not counts).
  S <- numeric(n_times)   # susceptible proportion
  I <- numeric(n_times)   # infected proportion
  R <- numeric(n_times)   # recovered/immune proportion

  # Initial conditions (day 1)
  S[1] <- 1 - I_init / pop   # everyone susceptible except seeded infections
  I[1] <- I_init / pop       # initial infected fraction
  R[1] <- 0                  # assume no one is recovered at start

  # --- Daily incidence (new infections per day, in COUNTS) ---
  incidence <- numeric(n_times)
  incidence[1] <- I[1] * pop   # first day: seeded infections

  # --- Main loop: update states day by day ---
  for (t in 2:n_times) {
    # Flows (proportion of population per day)
    inf_flow   <- beta_vec[t-1] * S[t-1] * I[t-1]   # new infections
    recov_flow <- gamma * I[t-1]                    # recoveries
    wane_flow  <- omega * R[t-1]                    # waning immunity

    # Update compartments
    S[t] <- S[t-1] - inf_flow + wane_flow
    I[t] <- I[t-1] + inf_flow - recov_flow
    R[t] <- R[t-1] + recov_flow - wane_flow

    # Prevent small numerical drift (e.g., negatives, sums not exactly 1)
    S[t] <- max(S[t], 0); I[t] <- max(I[t], 0); R[t] <- max(R[t], 0)
    tot <- S[t] + I[t] + R[t]
    if (!isTRUE(all.equal(tot, 1, tolerance = 1e-10))) {
      # renormalize so S+I+R = 1
      S[t] <- S[t] / tot
      I[t] <- I[t] / tot
      R[t] <- R[t] / tot
    }

    # Record daily incidence (counts, not proportions)
    incidence[t] <- pop * inf_flow
  }

  # --- Return results as a list ---
  list(
    time = seq_len(n_times),   # vector of days (1, 2, ..., n_times)
    S = S, I = I, R = R,       # proportions of population
    incidence = incidence,     # number of new infections per day
    params = list(             # record parameters for reproducibility
      n_times = n_times,
      pop     = pop,
      I_init  = I_init,
      gamma   = gamma,
      omega   = omega,
      beta    = beta_vec
    )
  )
}
