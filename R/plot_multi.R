# Deterministic: I(t) lines
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

# Deterministic: incidence lines
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

# Stochastic: I(t) ribbon + mean
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

# Stochastic: incidence ribbon + mean
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
