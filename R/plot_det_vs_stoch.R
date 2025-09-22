#' plot_det_vs_stoch
#' @title Overlay deterministic series over stochastic ribbon + mean
#' @description
#' Draws a quantile ribbon from many stochastic runs and overlays:
#'  - stochastic mean (black line)
#'  - deterministic series (colored line)
#'
#' @param det   Deterministic output from simulate_sirs_det().
#' @param stoch Stochastic output from simulate_sirs() (multi-run).
#' @param state One of "I" (proportion infected) or "incidence" (new cases).
#' @param probs Length-2 vector of quantiles for the ribbon, e.g., c(0.1, 0.9).
#' @param det_col Color for deterministic line.
#' @param main   Plot title (optional).
#'
#' @return Invisibly returns a list with mean and quantiles (for testing).
#' @examples
#' # plot_det_vs_stoch(det, st, state = "I")
plot_det_vs_stoch <- function(det, stoch,
                              state = c("I", "incidence"),
                              probs = c(0.1, 0.9),
                              det_col = "#0072B2",
                              main = NULL) {
  state <- match.arg(state)
  
  # --- basic checks ---
  if (is.null(det$time)) stop("det must have $time")
  if (is.null(stoch$proportions) || is.null(stoch$cases))
    stop("stoch must come from simulate_sirs() (needs $proportions and $cases).")
  
  # --- extract deterministic series and stochastic matrix [time x sims] ---
  x <- det$time
  if (state == "I") {
    if (is.null(det$I)) stop("det$I not found.")
    y_det <- det$I                                      # vector [time]
    y_mat <- stoch$proportions[, , "I", drop = FALSE]   # [time x sims x 1]
    y_mat <- y_mat[, , 1]                               # [time x sims]
    ylab <- "I (proportion)"
  } else {
    if (is.null(det$incidence)) stop("det$incidence not found.")
    y_det <- det$incidence
    y_mat <- stoch$cases                                # [time x sims]
    ylab <- "new infections (count)"
  }
  
  # --- compute band + mean across sims ---
  if (!is.numeric(probs) || length(probs) != 2) stop("probs must be length-2 numeric.")
  qfun <- function(M, probs) {
    qs <- apply(M, 1, quantile, probs = probs, na.rm = TRUE)
    list(low = qs[1, ], high = qs[2, ], mean = rowMeans(M, na.rm = TRUE))
  }
  band <- qfun(y_mat, probs)
  
  # --- helpers ---
  ribbon <- function(x, low, high, col = "gray", alpha = 0.25) {
    polygon(c(x, rev(x)), c(low, rev(high)), border = NA,
            col = grDevices::adjustcolor(col, alpha.f = alpha))
  }
  
  # --- y-limits to fit everything ---
  ylim <- range(y_det, band$low, band$high, 0, na.rm = TRUE)
  
  # --- draw ---
  plot(x, y_det, type = "n", xlab = "day", ylab = ylab, ylim = ylim,
       main = if (is.null(main)) sprintf("%s: det vs stochastic", state) else main)
  ribbon(x, band$low, band$high, col = "gray", alpha = 0.25)  # ribbon first
  lines(x, band$mean, lwd = 2, col = "black")                  # stochastic mean
  lines(x, y_det,    lwd = 2, col = det_col)                   # deterministic
  legend("topright",
         c(sprintf("stoch %d-%d%%", round(100*probs[1]), round(100*probs[2])),
           "stoch mean", "deterministic"),
         lty = 1, lwd = c(8, 2, 2),
         col = c(grDevices::adjustcolor("gray", 0.25), "black", det_col),
         bty = "n", seg.len = 2
  )
  
  invisible(list(mean = band$mean, low = band$low, high = band$high))
}
