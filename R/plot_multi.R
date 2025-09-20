#' @title Plot deterministic infected proportion by population
#' @description
#' Draw line plots of \eqn{I(t)} (infected proportion) for each population
#' from a deterministic multi-population simulation.
#'
#' @param sim A list-like simulation object containing:
#'   - `time` (numeric vector): time axis.
#'   - `I` (numeric matrix): infected proportions `[time x P]`.
#'
#' @return Invisibly returns `NULL`. Called for its side-effect of producing a plot.
#'
#' @examples
plot_multi_I <- function(sim) {
  stopifnot(!is.null(sim$I))
  Y <- sim$I
  P <- ncol(Y)
  
  matplot(sim$time, Y, type = "l", lwd = 2, lty = 1,
          col = seq_len(P),
          xlab = "day", ylab = "I (proportion)",
          main = "Infected by population")
  
  legend("topright", paste("Pop", seq_len(P)),
         col = seq_len(P), lty = 1, lwd = 2, bty = "n")
}

#' @title Plot deterministic incidence by population
#' @description
#' Draws line plots of daily incidence for each population from a deterministic
#' multi-population simulation. Optionally scales to cases per million.
#'
#' @param sim A list-like simulation object containing:
#'   - `time` (numeric vector): time axis.
#'   - `incidence` (numeric matrix): daily new infections `[time x P]`.
#'   - `params$pop_vec` (numeric vector, optional): population sizes for scaling.
#' @param per_million Logical. If `TRUE`, scale incidence to cases per million.
#'
#' @return
#' Invisibly returns `NULL`. Called for its side-effect of producing a plot.
#'
#' @examples
#' set.seed(2)
#' sim <- list(
#'   time = 1:60,
#'   incidence = matrix(rpois(60*3, lambda = 5), ncol = 3),
#'   params = list(pop_vec = c(1e6, 5e5, 2e6))
#' )
#' plot_multi_incidence(sim, per_million = TRUE)
#'
#' @importFrom graphics matplot legend
#' @export
plot_multi_incidence <- function(sim, per_million = FALSE) {
  stopifnot(!is.null(sim$incidence))
  Y <- sim$incidence
  P <- ncol(Y)
  
  if (per_million) {
    pops <- sim$params$pop_vec
    Y <- sweep(Y, 2, pops / 1e6, "/")
  }
  
  matplot(sim$time, Y, type = "l", lwd = 2, lty = 1,
          col = seq_len(P),
          xlab = "day",
          ylab = if (per_million) "cases per million" else "new infections",
          main = "Incidence by population")
  
  legend("topright", paste("Pop", seq_len(P)),
         col = seq_len(P), lty = 1, lwd = 2, bty = "n")
}

#' @title Plot stochastic infected proportion with uncertainty bands
#' @description
#' For each population, plots the mean \eqn{I(t)} (infected proportion) across
#' stochastic simulations with a ribbon showing central quantile bands.
#'
#' @param sim A list-like simulation object containing:
#'   - `time` (numeric vector): time axis.
#'   - `proportions` (numeric array): `[time x sims x P x state]`, including `"I"`.
#' @param probs Numeric length-2 vector of lower and upper quantiles (e.g., `c(0.05, 0.95)`).
#'
#' @return
#' Invisibly returns `NULL`. Called for its side-effect of producing plots.
#'
#' @examples
#' # Fake structure for demonstration:
#' set.seed(3)
#' time <- 1:50; sims <- 30; P <- 4
#' Iarr <- array(runif(length(time)*sims*P, 0, 0.05),
#'               dim = c(length(time), sims, P))
#' props <- array(NA_real_, dim = c(length(time), sims, P, 1),
#'                dimnames = list(NULL, NULL, NULL, c("I")))
#' props[,,, "I"] <- Iarr
#' sim <- list(time = time, proportions = props)
#' plot_multi_stoch_I(sim, probs = c(0.1, 0.9))
#'
#' @importFrom graphics par plot polygon lines legend
#' @importFrom grDevices adjustcolor
#' @importFrom stats quantile
#' @export
plot_multi_stoch_I <- function(sim, probs = c(0.05, 0.95)) {
  stopifnot(!is.null(sim$proportions))
  time <- sim$time
  P <- dim(sim$proportions)[3]
  
  # simple square-ish panel layout
  rows <- floor(sqrt(P)); cols <- ceiling(P / rows)
  op <- par(mfrow = c(rows, cols))
  on.exit(par(op), add = TRUE)
  
  # clamp probs to [0,1]
  probs <- pmax(pmin(probs, 1), 0)
  
  for (p in seq_len(P)) {
    I_mat <- sim$proportions[, , p, "I"]  # [time x sims]
    
    qL <- apply(I_mat, 1, quantile, probs = probs[1], na.rm = TRUE)
    qU <- apply(I_mat, 1, quantile, probs = probs[2], na.rm = TRUE)
    m  <- rowMeans(I_mat, na.rm = TRUE)
    
    plot(time, m, type = "n",
         ylim = range(c(qL, qU, m), na.rm = TRUE),
         xlab = "day", ylab = "I proportion",
         main = paste("Population", p))
    polygon(c(time, rev(time)), c(qL, rev(qU)),
            border = NA, col = adjustcolor("gray70", alpha.f = 0.6))
    lines(time, m, lwd = 2, col = "blue")
    
    legend("topright",
           c("Mean I(t)",
             sprintf("%d–%d%% band", round(probs[1]*100), round(probs[2]*100))),
           lty = c(1, NA), lwd = c(2, 10),
           col = c("blue", "gray70"), bty = "n")
  }
}

#' @title Plot stochastic incidence with uncertainty bands
#' @description
#' For each population, plots mean daily incidence across stochastic simulations
#' with a ribbon showing central quantile bands. Optionally scales to cases per million.
#'
#' @param sim A list-like simulation object containing:
#'   - `time` (numeric vector): time axis.
#'   - `cases` (numeric array): daily cases `[time x sims x P]`.
#'   - `params$pop_vec` (numeric vector, optional): population sizes for scaling.
#' @param per_million Logical. If `TRUE`, scale to cases per million.
#' @param probs Numeric length-2 vector of lower and upper quantiles.
#'
#' @return
#' Invisibly returns `NULL`. Called for its side-effect of producing plots.
#'
#' @examples
#' set.seed(4)
#' time <- 1:40; sims <- 20; P <- 3
#' cases <- array(rpois(length(time)*sims*P, 4), dim = c(length(time), sims, P))
#' sim <- list(time = time, cases = cases, params = list(pop_vec = c(8e5, 1.2e6, 2e6)))
#' plot_multi_stoch_incidence(sim, per_million = TRUE, probs = c(0.05, 0.95))
#'
#' @importFrom graphics par plot polygon lines legend
#' @importFrom grDevices adjustcolor
#' @importFrom stats quantile
#' @export
plot_multi_stoch_incidence <- function(sim, per_million = FALSE, probs = c(0.05, 0.95)) {
  stopifnot(!is.null(sim$cases))
  time <- sim$time
  P <- dim(sim$cases)[3]
  
  # simple square-ish panel layout
  rows <- floor(sqrt(P)); cols <- ceiling(P / rows)
  op <- par(mfrow = c(rows, cols))
  on.exit(par(op), add = TRUE)
  
  # clamp probs to [0,1]
  probs <- pmax(pmin(probs, 1), 0)
  
  scale_vec <- if (per_million) sim$params$pop_vec / 1e6 else rep(1, P)
  
  for (p in seq_len(P)) {
    Y <- sim$cases[, , p] / scale_vec[p]  # [time x sims]
    
    qL <- apply(Y, 1, quantile, probs = probs[1], na.rm = TRUE)
    qU <- apply(Y, 1, quantile, probs = probs[2], na.rm = TRUE)
    m  <- rowMeans(Y, na.rm = TRUE)
    
    plot(time, m, type = "n",
         ylim = range(c(qL, qU, m), na.rm = TRUE),
         xlab = "day",
         ylab = if (per_million) "cases per million" else "daily cases",
         main = paste("Population", p))
    polygon(c(time, rev(time)), c(qL, rev(qU)),
            border = NA, col = adjustcolor("gray70", alpha.f = 0.6))
    lines(time, m, lwd = 2, col = "blue")
    
    legend("topright",
           c("Mean cases",
             sprintf("%d–%d%% band", round(probs[1]*100), round(probs[2]*100))),
           lty = c(1, NA), lwd = c(2, 10),
           col = c("blue", "gray70"), bty = "n")
  }
}
