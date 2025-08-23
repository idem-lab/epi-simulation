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
