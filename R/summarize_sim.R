#' @title Summarize key outbreak metrics from a simulation
#' @description
#' Extract headline diagnostics from SIRS simulations: peak infected proportion
#' and its day, peak daily incidence (count) and its day, and the final recovered
#' proportion. Works for **deterministic single-run** outputs or **stochastic
#' single-population, multi-run** outputs.
#'
#' @param sim A simulation result list. See **Details** for the required fields.
#'
#' @return A `data.frame` with one row per simulation path and columns:
#' \describe{
#'   \item{sims}{Simulation index (1 for deterministic case).}
#'   \item{peak_I}{Maximum infected **proportion** over time.}
#'   \item{peak_I_day}{Time index (day) of `peak_I`.}
#'   \item{peak_incidence}{Maximum daily **count** of new infections.}
#'   \item{peak_incidence_day}{Time index (day) of `peak_incidence`.}
#'   \item{final_R}{Recovered/immune **proportion** at the last day.}
#' }
#'
#' @examples
#' # Deterministic example
#' # det <- simulate_sirs_det(n_times = 60, pop = 1e5, I_init = 10,
#' #                          beta = 0.16, gamma = 1/7, omega = 1/30)
#' # summarize_sim(det)
#'
#' # Stochastic (single population) example
#' # st <- simulate_sirs(n_times = 120, pop = 8e4, I_init = 12,
#' #                     beta = 0.16, gamma = 1/7, omega = 1/30,
#' #                     epsilon = 1e-4, n_sims = 25, seed = 42)
#' # head(summarize_sim(st))
#' 
#' @export
summarize_sim <- function(sim) {
  # This function extracts a few headline metrics from a simulation result:
  #  - peak_I:       maximum infected PROPORTION (or fraction) over time
  #  - peak_I_day:   day at which I(t) peaks
  #  - peak_incidence: maximum daily new infections (COUNT)
  #  - peak_incidence_day: day of that max incidence
  #  - final_R:      final recovered/immune PROPORTION at the last time step
  #
  # It supports two shapes of input:
  #  (A) Deterministic single-pop (from simulate_sirs_det): S/I/R vectors + incidence vector
  #  (B) Stochastic / multi-run (from simulate_sirs): proportions[time, sim, {S,I,R}] + cases[time, sim]

  # ---- Case A: deterministic single-run (vectors S/I/R in sim, no 'proportions' array) ----
  if (!is.null(sim$S) && is.null(sim$proportions)) {
    # Pull vectors
    S   <- sim$S
    I   <- sim$I
    R   <- sim$R
    inc <- sim$incidence

    # Basic guards against NA (just in case)
    I[!is.finite(I)]     <- NA_real_
    inc[!is.finite(inc)] <- NA_real_
    R[!is.finite(R)]     <- NA_real_

    # Peak infected (proportion) and the day it occurs
    peak_I   <- max(I, na.rm = TRUE)
    peak_day <- which.max(I)

    # Peak daily incidence (count) and the day it occurs
    peak_inc     <- max(inc, na.rm = TRUE)
    peak_inc_day <- which.max(inc)

    # Final recovered proportion (last day)
    final_R <- tail(R, 1)

    # Return a one-row data frame for consistency with the multi-sim case
    return(data.frame(
      sims = 1,
      peak_I = peak_I, peak_I_day = peak_day,
      peak_incidence = peak_inc, peak_incidence_day = peak_inc_day,
      final_R = final_R
    ))
  }

  # ---- Case B: stochastic / multi-run (array proportions[time, sim, {S,I,R}] + matrix cases[time, sim]) ----
  # Expect:
  #   sim$proportions: 3D array [time, sim, state], where state ∈ {"S","I","R"}
  #   sim$cases:       2D matrix [time, sim] of daily counts
  I_arr <- sim$proportions[, , "I", drop = FALSE]  # keep 3rd dim for clarity
  inc_m <- sim$cases
  n_sims <- dim(I_arr)[2]

  # Loop over sims and compute the same metrics for each
  out <- lapply(seq_len(n_sims), function(j) {
    Ii   <- I_arr[, j, 1]
    inci <- inc_m[, j]

    # Basic guards against NA (rare but safe)
    Ii[!is.finite(Ii)]     <- NA_real_
    inci[!is.finite(inci)] <- NA_real_

    data.frame(
      sims = j,
      peak_I = max(Ii, na.rm = TRUE),
      peak_I_day = which.max(Ii),
      peak_incidence = max(inci, na.rm = TRUE),
      peak_incidence_day = which.max(inci),
      final_R = sim$proportions[nrow(sim$proportions), j, "R"]
    )
  })

  # Bind rows into one data frame: one row per simulation
  do.call(rbind, out)
}
