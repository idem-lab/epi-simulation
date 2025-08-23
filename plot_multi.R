plot_multi_I <- function(sim) {
  # Plot the infected proportion (I) over time for each population in a multi-pop sim.
  # - sim: output from simulate_sirs_multi()
  # - sim$I: matrix [time x groups] of infected proportions

  stopifnot(!is.null(sim$I))

  matplot(sim$time, sim$I, type = "l", lwd = 2, lty = 1,
          col = seq_len(ncol(sim$I)),    # one color per population
          xlab = "day", ylab = "I (proportion)",
          main = "Infected by population")

  legend("topright", paste("Pop", seq_len(ncol(sim$I))),
         col = seq_len(ncol(sim$I)), lty = 1, lwd = 2, bty = "n")
}


plot_multi_incidence <- function(sim, per_million = FALSE) {
  # Plot daily incidence (new infections) for each population.
  # - sim: output from simulate_sirs_multi()
  # - sim$incidence: matrix [time x groups] of daily new infections (counts)
  # - per_million: if TRUE, normalize incidence per 1,000,000 people in each group

  stopifnot(!is.null(sim$incidence))

  Y <- sim$incidence

  # Optionally scale incidence per million population
  if (per_million) {
    pops <- sim$params$pop_vec
    Y <- sweep(Y, 2, pops / 1e6, "/")   # divide each column by (pop/1e6)
  }

  matplot(sim$time, Y, type = "l", lwd = 2, lty = 1,
          col = seq_len(ncol(Y)),       # one color per population
          xlab = "day",
          ylab = if (per_million) "cases per million" else "new infections",
          main = "Incidence by population")

  legend("topright", paste("Pop", seq_len(ncol(Y))),
         col = seq_len(ncol(Y)), lty = 1, lwd = 2, bty = "n")
}
