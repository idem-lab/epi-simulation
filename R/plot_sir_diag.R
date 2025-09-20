#' @title Plot SIR trajectories and/or incidence diagnostics
#' @description
#' Produce diagnostic plots from a single-population SIR/SIRS simulation:
#' S/I/R proportions, daily incidence, side-by-side panels, or an overlay with
#' a secondary axis for incidence.
#'
#' @param sim A list-like simulation object with numeric vectors:
#'   `time`, `S`, `I`, `R`, and `incidence` (all of equal length).
#' @param which Character. One of
#'   `c("both_side","overlay","sir","incidence")`. Selects the view to render.
#'
#' @return Invisibly returns `NULL`. Called for its side-effects of plotting.
#'
#' @examples
#' # Minimal reproducible example
#' set.seed(1)
#' t <- 1:120
#' S <- pmax(0, 1 - 0.6*(1 - exp(-t/60)))
#' I <- pmax(0, 0.15*exp(-(t-40)^2/800))
#' R <- pmax(0, 1 - S - I)
#' inc <- pmax(0, round(diff(c(0, I))*1000 + rnorm(length(t), 0, 2)))
#' sim <- list(time = t, S = S, I = I, R = R, incidence = inc)
#'
#' plot_sir_diag(sim, which = "sir")
#' plot_sir_diag(sim, which = "incidence")
#' plot_sir_diag(sim, which = "both_side")
#' plot_sir_diag(sim, which = "overlay")
#'
#' @importFrom graphics par plot lines legend axis mtext
#' @export
plot_sir_diag <- function(
    sim,
    which = c("both_side","overlay","sir","incidence")
) {
  # ---- choose which view to render ----
  which <- match.arg(which)

  # ---- minimal input validation (require the usual fields) ----
  need <- c("time","S","I","R","incidence")
  miss <- setdiff(need, names(sim))
  if (length(miss)) stop("sim is missing: ", paste(miss, collapse=", "))

  # ---------------------------
  # 1) Only S, I, R panel
  # ---------------------------
  if (which == "sir") {
    # Plot S, I, R as proportions on the same axes (y in [0,1])
    plot(sim$time, sim$S, type="l", lwd=2, ylim=c(0,1), col="#0072B2",
         xlab="day", ylab="proportion", main="S, I, R")
    lines(sim$time, sim$I, lwd=2, col="#CC79A7")
    lines(sim$time, sim$R, lwd=2, col="#009E73")
    legend("topright", c("S","I","R"),
           lty=1, lwd=2, col=c("#0072B2","#CC79A7","#009E73"), bty="n")
    return(invisible(NULL))
  }

  # ---------------------------
  # 2) Only incidence panel
  # ---------------------------
  if (which == "incidence") {
    # Plot daily new infections (counts)
    plot(sim$time, sim$incidence, type="l", lwd=2, col="darkgray",
         xlab="day", ylab="new infections (count)", main="Daily incidence")
    return(invisible(NULL))
  }

  # ---------------------------
  # 3) Side-by-side panels
  # ---------------------------
  if (which == "both_side") {
    # Two plots next to each other; restore par() afterwards
    oldpar <- par(mfrow = c(1,2))
    on.exit(par(oldpar), add = TRUE)

    # left: S, I, R (proportions)
    plot(sim$time, sim$S, type="l", lwd=2, ylim=c(0,1), col="#0072B2",
         xlab="day", ylab="proportion", main="S, I, R")
    lines(sim$time, sim$I, lwd=2, col="#CC79A7")
    lines(sim$time, sim$R, lwd=2, col="#009E73")
    legend("topright", c("S","I","R"),
           lty=1, lwd=2, col=c("#0072B2","#CC79A7","#009E73"), bty="n")

    # right: incidence (counts)
    plot(sim$time, sim$incidence, type="l", lwd=2, col="darkgray",
         xlab="day", ylab="new infections (count)", main="Daily incidence")
    return(invisible(NULL))
  }

  # ---------------------------
  # 4) Overlay incidence on S/I/R with a 2nd y-axis
  # ---------------------------
  if (which == "overlay") {
    # Base plot: S, I, R on [0,1] left y-axis
    plot(sim$time, sim$S, type="l", lwd=2, ylim=c(0,1), col="#0072B2",
         xlab="day", ylab="proportion", main="S, I, R + incidence")
    lines(sim$time, sim$I, lwd=2, col="#CC79A7")
    lines(sim$time, sim$R, lwd=2, col="#009E73")
    legend("topleft", c("S","I","R"),
           lty=1, lwd=2, col=c("#0072B2","#CC79A7","#009E73"), bty="n")

    # Add incidence as a second (right) axis.
    # par(new=TRUE) tells R to draw on top of the existing plot area.
    on.exit(par(new = FALSE), add = TRUE)
    par(new = TRUE)

    # Use the data's own range for the right axis limits
    rng <- range(sim$incidence, na.rm = TRUE)

    # Draw incidence without axes; we’ll add a right-side axis next
    plot(sim$time, sim$incidence, type="l", lwd=2, axes=FALSE, xlab="", ylab="",
         col="darkgray", ylim = rng)

    # Right-side axis + label
    axis(side = 4)
    mtext("incidence (count)", side = 4, line = 3)
    return(invisible(NULL))
  }
}
