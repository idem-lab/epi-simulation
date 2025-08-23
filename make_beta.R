make_beta <- function(
    n_times,                        # total simulation length (days)
    mode = c("constant","seasonal"),# choose type of beta: constant OR seasonal
    value = 0.16,                   # constant beta value (if mode = "constant")
    base = 0.16, amplitude = 0.2, phase = 0  # seasonal parameters
) {
  # Pick the mode (only one valid choice)
  mode <- match.arg(mode)

  if (mode == "constant") {
    # Same beta every day
    return(rep(value, n_times))
  } else {
    # Seasonal: beta oscillates like a cosine curve over a 365-day period
    # Formula: beta_t = base * (1 + amplitude * cos(2*pi*(t + phase)/365))
    # - base: average transmission rate
    # - amplitude: how much it fluctuates (0 = flat, closer to 1 = bigger swings)
    # - phase: shifts the seasonal curve (e.g., peak in winter vs summer)
    if (amplitude < 0 || amplitude >= 1) stop("amplitude must be in [0,1).")
    t <- seq_len(n_times)
    b <- base * (1 + amplitude * cos(2*pi*(t + phase)/365))

    # Just in case: prevent negative beta values
    b[b < 0] <- 0
    return(b)
  }
}
