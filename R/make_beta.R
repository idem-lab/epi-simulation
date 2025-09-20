#' @title Construct a time-varying transmission rate beta_t
#' 
#' @description
#' Returns a length-`n_times` vector of daily transmission rates \eqn{\beta_t}
#' for SIR/SIRS-style simulations. Supports a constant rate or a seasonal
#' cosine pattern over a 365-day period.
#'
#' @param n_times Integer. Total simulation length (days).
#' @param mode Character. One of `"constant"` or `"seasonal"`.
#' @param value Numeric. Constant \eqn{\beta} when `mode = "constant"`.
#' @param base Numeric. Average level when `mode = "seasonal"`.
#' @param amplitude Numeric in \eqn{[0,1)}. Seasonal fluctuation size.
#' @param phase Numeric. Phase shift (days) for seasonality.
#'
#' @return Numeric vector of length `n_times` giving \eqn{\beta_t}.
#'
#' @examples
#' # Constant beta
#' make_beta(n_times = 30, mode = "constant", value = 0.16)
#'
#' # Seasonal beta
#' make_beta(n_times = 365, mode = "seasonal", base = 0.18, amplitude = 0.25, phase = 30)
#'
#' @export
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
