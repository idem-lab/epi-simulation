# =========================================
# 6 plot version
# =========================================
plot_dashboard <- function(sim, probs = c(0.05, 0.95), per_million = FALSE,
                                    main = "Multi-pop dashboard") {
  stopifnot(!is.null(sim$time))
  # helpers
  qfun <- function(x, p) stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  safe_lab <- function(x) if (is.null(x)) "—" else if (length(x) > 6) paste0("[", length(x), " values]") else paste(x, collapse = ", ")
  has_det <- !is.null(sim$I) && is.matrix(sim$I)
  has_sto <- !is.null(sim$proportions) && length(dim(sim$proportions)) == 4
  P <- if (has_det) ncol(sim$I) else dim(sim$proportions)[3]
  time <- sim$time
  pops <- if (!is.null(sim$params$pop_vec)) sim$params$pop_vec else rep(1, P)
  
  # build series for I(t) and incidence overviews
  if (has_det) {
    I_over <- sim$I                                   # [time x P]
    Inc_over <- sim$incidence                         # [time x P]
  } else if (has_sto) {
    # mean across sims
    I_over   <- apply(sim$proportions[, , , "I"], c(1, 3), mean, na.rm = TRUE)  # [time x P]
    Inc_over <- apply(sim$cases, c(1, 3), mean, na.rm = TRUE)                    # [time x P]
  } else stop("Unsupported sim shape for multi-pop dashboard.")
  
  if (per_million) {
    I_ylab <- "I proportion (per million population)"
    Inc_ylab <- "daily cases (per million)"
    I_over   <- sweep(I_over,   2, pops / 1e6, "/")
    Inc_over <- sweep(Inc_over, 2, pops / 1e6, "/")
  } else {
    I_ylab <- "I proportion"
    Inc_ylab <- "daily cases"
  }
  
  # β handling
  beta_obj <- sim$params$beta
  beta_is_vec <- is.vector(beta_obj) && length(beta_obj) == length(time)
  beta_is_mat <- is.matrix(beta_obj) && nrow(beta_obj) == length(time) && ncol(beta_obj) == P
  
  # layout 3x2
  op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 2) + 0.1)
  on.exit(par(op), add = TRUE)
  
  # (1) I(t) overview
  matplot(time, I_over, type = "l", lwd = 2, lty = 1, col = seq_len(P),
          xlab = "day", ylab = I_ylab, main = paste(main, "- I(t) overview"))
  legend("topright", paste("Pop", seq_len(P)), col = seq_len(P), lty = 1, lwd = 2, bty = "n")
  
  # (2) Incidence overview
  matplot(time, Inc_over, type = "l", lwd = 2, lty = 1, col = seq_len(P),
          xlab = "day", ylab = Inc_ylab, main = "Incidence overview")
  legend("topright", paste("Pop", seq_len(P)), col = seq_len(P), lty = 1, lwd = 2, bty = "n")
  
  # (3) β(t)
  if (beta_is_vec) {
    plot(time, beta_obj, type = "l", lwd = 2, xlab = "day", ylab = expression(beta),
         main = expression(beta~"(t)"))
  } else if (beta_is_mat) {
    matplot(time, beta_obj, type = "l", lwd = 2, lty = 1, col = seq_len(P),
            xlab = "day", ylab = expression(beta), main = expression(beta~"(t) by population"))
    legend("topright", paste("Pop", seq_len(P)), col = seq_len(P), lty = 1, lwd = 2, bty = "n")
  } else {
    plot.new(); title(expression(beta~"(t) unavailable")); box()
  }
  
  # (4) Contact matrix heatmap (if available)
  C <- sim$params$C
  if (is.matrix(C) && nrow(C) == P && ncol(C) == P) {
    # pretty heatmap using image()
    image(z = t(C[nrow(C):1, ]), axes = FALSE, main = "Contact matrix (rows recv, cols src)")
    axis(1, at = seq(0, 1, length.out = P), labels = seq_len(P))
    axis(2, at = seq(0, 1, length.out = P), labels = rev(seq_len(P)))
    box()
  } else {
    plot.new(); title("No contact matrix provided"); box()
  }
  
  # (5) Peak I by population (with stochastic bands if available)
  if (has_det) {
    peak_I <- apply(sim$I, 2, max, na.rm = TRUE)
    bp <- barplot(peak_I, col = "gray80", border = "gray40",
                  main = "Peak I by population", ylab = "peak I (proportion)")
    text(bp, peak_I, labels = sprintf("%.3f", peak_I), pos = 3, cex = 0.9)
  } else {
    # stochastic: compute per-sim peaks then summarise by quantiles per pop
    I_arr <- sim$proportions[, , , "I"]     # [time x sims x P]
    n_sims <- dim(I_arr)[2]
    peak_mat <- apply(I_arr, c(2, 3), max, na.rm = TRUE)  # [sims x P]
    mean_p <- colMeans(peak_mat, na.rm = TRUE)
    qL <- apply(peak_mat, 2, qfun, p = probs[1])
    qU <- apply(peak_mat, 2, qfun, p = probs[2])
    bp <- barplot(mean_p, col = "gray80", border = "gray40",
                  main = sprintf("Peak I by population (n=%d sims)", n_sims),
                  ylab = sprintf("mean peak I (%.0f–%.0f%% band)", probs[1]*100, probs[2]*100))
    # error bars
    arrows(x0 = bp, y0 = qL, x1 = bp, y1 = qU, angle = 90, code = 3, length = 0.05, lwd = 2)
    text(bp, mean_p, labels = sprintf("%.3f", mean_p), pos = 3, cex = 0.9)
  }
  
  # (6) Parameter box
  plot.new(); box(); title("Parameters", line = 0.5)
  usr <- par("usr"); xmid <- (usr[1] + usr[2]) / 2
  y <- 0.92
  add_line <- function(lbl, val) { text(xmid, y, paste0(lbl, ": ", val), cex = 0.9); y <<- y - 0.095 }
  p <- sim$params
  add_line("n_times", length(time))
  add_line("P", P)
  add_line("pop_vec", safe_lab(p$pop_vec))
  add_line("I_init",  safe_lab(p$I_init))
  add_line("gamma",   if (!is.null(p$gamma)) signif(p$gamma, 4) else "—")
  add_line("omega",   if (!is.null(p$omega)) signif(p$omega, 4) else "—")
  add_line("epsilon", if (!is.null(p$epsilon)) signif(p$epsilon, 4) else "—")
  add_line("alpha",   if (!is.null(p$alpha))   signif(p$alpha, 4) else "—")
  add_line("n_sims",  if (!is.null(p$n_sims))  p$n_sims else (if (has_det) 1 else "—"))
  add_line("stochastic", if (!is.null(p$stochastic)) p$stochastic else has_sto)
}


# =========================================
# 9 plots version
# =========================================
plot_dashboard_v2 <- function(sim, probs = c(0.05, 0.95), per_million = FALSE,
                                       main = "Multi-pop dashboard") {
  stopifnot(!is.null(sim$time))
  qfun <- function(x, p) stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  
  has_det <- !is.null(sim$I) && is.matrix(sim$I)                    # det: S/I/R matrices [time x P]
  has_sto <- !is.null(sim$proportions) && length(dim(sim$proportions)) == 4   # stoch: props[time, sims, P, state]
  P <- if (has_det) ncol(sim$I) else dim(sim$proportions)[3]
  time <- sim$time
  pops <- if (!is.null(sim$params$pop_vec)) sim$params$pop_vec else rep(1, P)
  
  # helpers to produce means & ribbons across sims
  mean_state <- function(state) apply(sim$proportions[, , , state], c(1, 3), mean, na.rm = TRUE)
  qL_state   <- function(state) apply(sim$proportions[, , , state], c(1, 3), qfun, p = probs[1])
  qU_state   <- function(state) apply(sim$proportions[, , , state], c(1, 3), qfun, p = probs[2])
  
  # Build S/I/R time series (det → directly; stoch → mean)
  if (has_det) {
    S_over <- sim$S; I_over <- sim$I; R_over <- sim$R                     # [time x P]
    Inc_over <- sim$incidence                                             # [time x P]
  } else if (has_sto) {
    S_over <- mean_state("S"); I_over <- mean_state("I"); R_over <- mean_state("R")
    Inc_over <- apply(sim$cases, c(1, 3), mean, na.rm = TRUE)             # [time x P]
  } else stop("Unsupported sim shape.")
  
  # Optional per-million scaling for incidence
  Inc_ylab <- if (per_million) "daily cases (per million)" else "daily cases"
  if (per_million) Inc_over <- sweep(Inc_over, 2, pops / 1e6, "/")
  
  # Beta handling
  beta_obj <- sim$params$beta
  beta_is_vec <- is.vector(beta_obj) && length(beta_obj) == length(time)
  beta_is_mat <- is.matrix(beta_obj) && nrow(beta_obj) == length(time) && ncol(beta_obj) == P
  
  # Ribbon bounds
  if (has_sto) {
    I_qL <- qL_state("I"); I_qU <- qU_state("I")                         # [time x P]
    Inc_qL <- apply(sim$cases, c(1, 3), qfun, p = probs[1])
    Inc_qU <- apply(sim$cases, c(1, 3), qfun, p = probs[2])
    if (per_million) { Inc_qL <- sweep(Inc_qL, 2, pops / 1e6, "/"); Inc_qU <- sweep(Inc_qU, 2, pops / 1e6, "/") }
  }
  
  cols <- seq_len(P)                        
  alpha <- function(col, a = 0.35) grDevices::adjustcolor(col, alpha.f = a)
  
  # 3 x 3 [S, I(+ribbons), R] / [Inc(+ribbons), Beta, Contact] / [PeakI, PeakDay, Params]
  op <- par(mfrow = c(3, 3), mar = c(4, 4, 3, 2) + 0.1)
  on.exit(par(op), add = TRUE)
  
  # (1) S(t) overview
  matplot(time, S_over, type = "l", lwd = 2, lty = 1, col = cols,
          xlab = "day", ylab = "S proportion", main = paste(main, "- S(t) overview"))
  legend("topright", paste("Pop", seq_len(P)), col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (2) I(t) overview + ribbons (stoch only)
  ylim_I <- range(I_over, if (has_sto) c(I_qL, I_qU) else NULL, na.rm = TRUE)
  plot(time, I_over[,1], type = "n", ylim = ylim_I, xlab = "day", ylab = "I proportion",
       main = "I(t) overview")
  if (has_sto) {
    for (p in seq_len(P)) {
      polygon(c(time, rev(time)), c(I_qL[,p], rev(I_qU[,p])), border = NA, col = alpha(cols[p], 0.25))
    }
  }
  for (p in seq_len(P)) lines(time, I_over[,p], lwd = 2, col = cols[p])
  legend("topright", paste("Pop", seq_len(P)), col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (3) R(t) overview
  matplot(time, R_over, type = "l", lwd = 2, lty = 1, col = cols,
          xlab = "day", ylab = "R proportion", main = "R(t) overview")
  
  # (4) Incidence overview + ribbons
  ylim_inc <- range(Inc_over, if (has_sto) c(Inc_qL, Inc_qU) else NULL, na.rm = TRUE)
  plot(time, Inc_over[,1], type = "n", ylim = ylim_inc, xlab = "day", ylab = Inc_ylab,
       main = "Incidence overview")
  if (has_sto) {
    for (p in seq_len(P)) {
      polygon(c(time, rev(time)), c(Inc_qL[,p], rev(Inc_qU[,p])), border = NA, col = alpha(cols[p], 0.25))
    }
  }
  for (p in seq_len(P)) lines(time, Inc_over[,p], lwd = 2, col = cols[p])
  legend("topright", paste("Pop", seq_len(P)), col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (5) Beta(t)
  if (beta_is_vec) {
    plot(time, beta_obj, type = "l", lwd = 2, xlab = "day", ylab = expression(beta), main = expression(beta~"(t)"))
  } else if (beta_is_mat) {
    matplot(time, beta_obj, type = "l", lwd = 2, lty = 1, col = cols,
            xlab = "day", ylab = expression(beta), main = expression(beta~"(t) by population"))
    legend("topright", paste("Pop", seq_len(P)), col = cols, lty = 1, lwd = 2, bty = "n")
  } else { plot.new(); title(expression(beta~"(t) unavailable")); box() }
  
  # (6) Contact matrix
  C <- sim$params$C
  if (is.matrix(C) && nrow(C) == P && ncol(C) == P) {
    image(z = t(C[nrow(C):1, ]), axes = FALSE, main = "Contact matrix (rows recv, cols src)")
    axis(1, at = seq(0, 1, length.out = P), labels = seq_len(P))
    axis(2, at = seq(0, 1, length.out = P), labels = rev(seq_len(P)))
    box()
  } else { plot.new(); title("No contact matrix"); box() }
  
  # (7) Peak I by population（mean + band）
  if (has_det) {
    peak_I <- apply(sim$I, 2, max, na.rm = TRUE)
    bp <- barplot(peak_I, col = "gray80", border = "gray40", ylim = c(0, max(peak_I)*1.2),
                  main = "Peak I by population", ylab = "peak I (proportion)")
    text(bp, peak_I, labels = sprintf("%.3f", peak_I), pos = 3, cex = 0.9)
  } else {
    I_arr <- sim$proportions[, , , "I"]; n_s <- dim(I_arr)[2]
    peak_mat <- apply(I_arr, c(2, 3), max, na.rm = TRUE)         # [sims x P]
    mean_p <- colMeans(peak_mat, na.rm = TRUE)
    qL_p <- apply(peak_mat, 2, qfun, p = probs[1])
    qU_p <- apply(peak_mat, 2, qfun, p = probs[2])
    bp <- barplot(mean_p, col = "gray80", border = "gray40",
                  main = sprintf("Peak I by population (n=%d sims)", n_s),
                  ylab = sprintf("mean peak I (%.0f–%.0f%% band)", probs[1]*100, probs[2]*100),
                  ylim = c(0, max(qU_p)*1.15))
    arrows(bp, qL_p, bp, qU_p, angle = 90, code = 3, length = 0.05, lwd = 2)
    text(bp, mean_p, labels = sprintf("%.3f", mean_p), pos = 3, cex = 0.9)
  }
  
  # (8) Stochasticity magnitude
  if (has_sto) {
    width_I <- colMeans(I_qU - I_qL, na.rm = TRUE)          
    barplot(width_I, col = alpha("gray40", 0.6), border = "gray20",
            main = "Stochasticity (avg ribbon width of I)", ylab = "avg (qU - qL)")
  } else {
    plot.new(); title("Stochasticity (det only: none)"); box()
  }
  
  # (9) Parameters
  plot.new(); box(); title("Parameters", line = 0.5)
  usr <- par("usr"); xmid <- (usr[1] + usr[2]) / 2; y <- 0.92
  add_line <- function(lbl, val) { text(xmid, y, paste0(lbl, ": ", val), cex = 0.9); y <<- y - 0.095 }
  p <- sim$params
  add_line("n_times", length(time))
  add_line("P", P)
  add_line("pop_vec", if (is.null(p$pop_vec)) "—" else paste(p$pop_vec, collapse = ", "))
  add_line("I_init",  if (is.null(p$I_init))  "—" else paste(p$I_init, collapse = ", "))
  add_line("beta",    if (is.null(p$beta))    "—" else if (is.matrix(p$beta)) "[matrix]" else if (length(p$beta)>6) "[vector]" else paste(p$beta, collapse=","))
  add_line("gamma",   if (is.null(p$gamma))   "—" else signif(p$gamma, 4))
  add_line("1/gamma (days)", if (is.null(p$gamma)) "—" else sprintf("%.2f", 1/p$gamma))
  add_line("omega",   if (is.null(p$omega))   "—" else signif(p$omega, 4))
  add_line("1/omega (days)", if (is.null(p$omega)) "—" else sprintf("%.2f", 1/p$omega))
  add_line("epsilon", if (is.null(p$epsilon)) "—" else signif(p$epsilon, 4))
  add_line("alpha",   if (is.null(p$alpha))   "—" else signif(p$alpha, 4))
  add_line("n_sims",  if (!is.null(p$n_sims)) p$n_sims else (if (has_det) 1 else "—"))
  add_line("stochastic", if (!is.null(p$stochastic)) p$stochastic else has_sto)
}
