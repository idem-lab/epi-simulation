# ============================================================
# plot_dashboard.R — simple 6-plot & 9-plot dashboards
# (no contact panel switching; fixed colours; robust y-limits)
# ============================================================

# ---------- shared helpers ----------
.qfun <- function(x, p) stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE)
.safe_lab <- function(x) if (is.null(x)) "—" else if (length(x) > 6) paste0("[", length(x), " values]") else paste(x, collapse = ", ")
.safe_range <- function(...) {
  r <- range(..., na.rm = TRUE)
  if (!all(is.finite(r))) return(c(0, 1))
  if (r[1] == r[2]) {
    bump <- ifelse(r[2] == 0, 1, abs(r[2]) * 0.05)
    return(r + c(-bump, bump))
  }
  r
}
# EXACT palette requested
.pop_cols <- function(P) {
  base <- c("#000000", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#56B4E9", "#F0E442")
  base[seq_len(min(P, length(base)))]
}

# ============================================================
# 6-plot dashboard (overlay in panel 2)
# Panels: I, DetVsStoch, Beta, Incidence, PeakI, Params
# ============================================================
plot_dashboard <- function(sim, probs = c(0.05, 0.95), per_million = FALSE,
                           main = "Multi-pop dashboard",
                           det_stoch = NULL) {
  stopifnot(!is.null(sim$time))
  
  has_det <- !is.null(sim$I) && is.matrix(sim$I)
  has_sto <- !is.null(sim$proportions) && length(dim(sim$proportions)) == 4
  P <- if (has_det) ncol(sim$I) else dim(sim$proportions)[3]
  time <- sim$time
  pops <- if (!is.null(sim$params$pop_vec)) sim$params$pop_vec else rep(1, P)
  pop_labels <- if (!is.null(sim$params$pop_names)) sim$params$pop_names else paste("Pop", seq_len(P))
  cols <- .pop_cols(P)
  
  # series
  if (has_det) {
    I_over <- sim$I
    Inc_over <- sim$incidence
  } else if (has_sto) {
    I_over   <- apply(sim$proportions[, , , "I"], c(1, 3), mean, na.rm = TRUE)
    Inc_over <- apply(sim$cases, c(1, 3), mean, na.rm = TRUE)
  } else stop("Unsupported sim shape.")
  
  if (per_million) {
    I_ylab <- "I proportion (per million population)"
    Inc_ylab <- "daily cases (per million)"
    I_over   <- sweep(I_over,   2, pops / 1e6, "/")
    Inc_over <- sweep(Inc_over, 2, pops / 1e6, "/")
  } else { I_ylab <- "I proportion"; Inc_ylab <- "daily cases" }
  
  beta_obj <- sim$params$beta
  beta_is_vec <- is.vector(beta_obj) && length(beta_obj) == length(time)
  beta_is_mat <- is.matrix(beta_obj) && nrow(beta_obj) == length(time) && ncol(beta_obj) == P
  
  op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 2) + 0.1); on.exit(par(op), add = TRUE)
  
  # (1) I(t)
  matplot(time, I_over, type = "l", lwd = 2, lty = 1, col = cols,
          xlab = "day", ylab = I_ylab, main = paste(main, " - I(t) overview"))
  legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (2) Det vs Stoch overlay (or fallback to incidence)
  if (!is.null(det_stoch) && exists("plot_det_vs_stoch", mode = "function")) {
    do.call(plot_det_vs_stoch, args = c(det_stoch, list()))
  } else {
    matplot(time, Inc_over, type = "l", lwd = 2, lty = 1, col = cols,
            xlab = "day", ylab = Inc_ylab, main = "Incidence overview",
            ylim = .safe_range(Inc_over))
    legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  }
  
  # (3) β(t)
  if (beta_is_vec) {
    plot(time, beta_obj, type = "l", lwd = 2, xlab = "day", ylab = expression(beta), main = expression(beta~"(t)"))
  } else if (beta_is_mat) {
    matplot(time, beta_obj, type = "l", lwd = 2, lty = 1, col = cols,
            xlab = "day", ylab = expression(beta), main = expression(beta~"(t) by population"))
    legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  } else { plot.new(); title(expression(beta~"(t) unavailable")); box() }
  
  # (4) Incidence overview (always available)
  matplot(time, Inc_over, type = "l", lwd = 2, lty = 1, col = cols,
          xlab = "day", ylab = Inc_ylab, main = "Incidence overview",
          ylim = .safe_range(Inc_over))
  legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (5) Peak I
  if (has_det) {
    peak_I <- apply(sim$I, 2, max, na.rm = TRUE)
    bp <- barplot(peak_I, col = "gray80", border = "gray40", ylim = .safe_range(peak_I),
                  main = "Peak I by population", ylab = "peak I (proportion)")
    text(bp, peak_I, labels = sprintf("%.3f", peak_I), pos = 3, cex = 0.9)
  } else {
    I_arr <- sim$proportions[, , , "I"]; n_sims <- dim(I_arr)[2]
    peak_mat <- apply(I_arr, c(2, 3), max, na.rm = TRUE)
    mean_p <- colMeans(peak_mat, na.rm = TRUE)
    qL <- apply(peak_mat, 2, .qfun, p = probs[1])
    qU <- apply(peak_mat, 2, .qfun, p = probs[2])
    bp <- barplot(mean_p, col = "gray80", border = "gray40",
                  main = sprintf("Peak I by population (n=%d sims)", n_sims),
                  ylab = sprintf("mean peak I (%.0f–%.0f%% band)", probs[1]*100, probs[2]*100),
                  ylim = .safe_range(0, qU))
    arrows(bp, qL, bp, qU, angle = 90, code = 3, length = 0.05, lwd = 2)
    text(bp, mean_p, labels = sprintf("%.3f", mean_p), pos = 3, cex = 0.9)
  }
  
  # (6) Parameters
  plot.new(); box(); title("Parameters", line = 0.5)
  usr <- par("usr"); xmid <- (usr[1] + usr[2]) / 2; y <- 0.92
  add_line <- function(lbl, val) { text(xmid, y, paste0(lbl, ": ", val), cex = 0.9); y <<- y - 0.095 }
  p <- sim$params
  add_line("n_times", length(time)); add_line("P", P)
  add_line("pop_vec", .safe_lab(p$pop_vec)); add_line("I_init", .safe_lab(p$I_init))
  add_line("beta",    if (is.null(p$beta)) "—" else if (is.matrix(p$beta)) "[matrix]" else if (length(p$beta)>6) "[vector]" else paste(p$beta, collapse=","))
  add_line("gamma",   if (!is.null(p$gamma)) signif(p$gamma, 4) else "—")
  add_line("1/gamma (days)", if (!is.null(p$gamma)) sprintf("%.2f", 1/p$gamma) else "—")
  add_line("omega",   if (!is.null(p$omega)) signif(p$omega, 4) else "—")
  add_line("1/omega (days)", if (!is.null(p$omega)) sprintf("%.2f", 1/p$omega) else "—")
  add_line("epsilon", if (!is.null(p$epsilon)) signif(p$epsilon, 4) else "—")
  add_line("alpha",   if (!is.null(p$alpha))   signif(p$alpha, 4) else "—")
  add_line("n_sims",  if (!is.null(p$n_sims))  p$n_sims else (if (has_det) 1 else "—"))
  add_line("stochastic", if (!is.null(p$stochastic)) p$stochastic else has_sto)
}

# ============================================================
# 9-plot dashboard (overlay in panel 6 — middle-right)
# Panels: S, I(+ribbons), R, Inc(+ribbons), Beta, DetVsStoch, PeakI, StochMag, Params
# ============================================================
plot_dashboard_v2 <- function(sim, probs = c(0.05, 0.95), per_million = FALSE,
                              main = "Multi-pop dashboard",
                              det_stoch = NULL) {
  stopifnot(!is.null(sim$time))
  
  has_det <- !is.null(sim$I) && is.matrix(sim$I)
  has_sto <- !is.null(sim$proportions) && length(dim(sim$proportions)) == 4
  P <- if (has_det) ncol(sim$I) else dim(sim$proportions)[3]
  time <- sim$time
  pops <- if (!is.null(sim$params$pop_vec)) sim$params$pop_vec else rep(1, P)
  pop_labels <- if (!is.null(sim$params$pop_names)) sim$params$pop_names else paste("Pop", seq_len(P))
  cols <- .pop_cols(P)
  alpha_col <- function(col, a = 0.35) grDevices::adjustcolor(col, alpha.f = a)
  
  # series
  if (has_det) {
    S_over <- sim$S; I_over <- sim$I; R_over <- sim$R
    Inc_over <- sim$incidence
  } else if (has_sto) {
    S_over <- apply(sim$proportions[, , , "S"], c(1, 3), mean, na.rm = TRUE)
    I_over <- apply(sim$proportions[, , , "I"], c(1, 3), mean, na.rm = TRUE)
    R_over <- apply(sim$proportions[, , , "R"], c(1, 3), mean, na.rm = TRUE)
    Inc_over <- apply(sim$cases, c(1, 3), mean, na.rm = TRUE)
  } else stop("Unsupported sim shape.")
  
  Inc_ylab <- if (per_million) "daily cases (per million)" else "daily cases"
  if (per_million) Inc_over <- sweep(Inc_over, 2, pops / 1e6, "/")
  
  beta_obj <- sim$params$beta
  beta_is_vec <- is.vector(beta_obj) && length(beta_obj) == length(time)
  beta_is_mat <- is.matrix(beta_obj) && nrow(beta_obj) == length(time) && ncol(beta_obj) == P
  
  if (has_sto) {
    I_qL <- apply(sim$proportions[, , , "I"], c(1, 3), .qfun, p = probs[1])
    I_qU <- apply(sim$proportions[, , , "I"], c(1, 3), .qfun, p = probs[2])
    Inc_qL <- apply(sim$cases, c(1, 3), .qfun, p = probs[1])
    Inc_qU <- apply(sim$cases, c(1, 3), .qfun, p = probs[2])
    if (per_million) { Inc_qL <- sweep(Inc_qL, 2, pops / 1e6, "/"); Inc_qU <- sweep(Inc_qU, 2, pops / 1e6, "/") }
  }
  
  op <- par(mfrow = c(3, 3), mar = c(4, 4, 3, 2) + 0.1); on.exit(par(op), add = TRUE)
  
  # (1) S(t)
  matplot(time, S_over, type = "l", lwd = 2, lty = 1, col = cols,
          xlab = "day", ylab = "S proportion", main = paste(main, " - S(t) overview"))
  legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (2) I(t)
  ylim_I <- if (has_sto) .safe_range(I_over, I_qL, I_qU) else .safe_range(I_over)
  plot(time, I_over[,1], type = "n", ylim = ylim_I, xlab = "day", ylab = "I proportion",
       main = "I(t) overview")
  if (has_sto) for (p in seq_len(P)) polygon(c(time, rev(time)), c(I_qL[,p], rev(I_qU[,p])), border = NA, col = alpha_col(cols[p], 0.25))
  for (p in seq_len(P)) lines(time, I_over[,p], lwd = 2, col = cols[p])
  legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (3) R(t)
  matplot(time, R_over, type = "l", lwd = 2, lty = 1, col = cols,
          xlab = "day", ylab = "R proportion", main = "R(t) overview")
  
  # (4) Incidence
  ylim_inc <- if (has_sto) .safe_range(Inc_over, Inc_qL, Inc_qU) else .safe_range(Inc_over)
  plot(time, Inc_over[,1], type = "n", ylim = ylim_inc, xlab = "day", ylab = Inc_ylab,
       main = "Incidence overview")
  if (has_sto) for (p in seq_len(P)) polygon(c(time, rev(time)), c(Inc_qL[,p], rev(Inc_qU[,p])), border = NA, col = alpha_col(cols[p], 0.25))
  for (p in seq_len(P)) lines(time, Inc_over[,p], lwd = 2, col = cols[p])
  legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  
  # (5) β(t)
  if (beta_is_vec) {
    plot(time, beta_obj, type = "l", lwd = 2, xlab = "day", ylab = expression(beta), main = expression(beta~"(t)"))
  } else if (beta_is_mat) {
    matplot(time, beta_obj, type = "l", lwd = 2, lty = 1, col = cols,
            xlab = "day", ylab = expression(beta), main = expression(beta~"(t) by population"))
    legend("topright", pop_labels, col = cols, lty = 1, lwd = 2, bty = "n")
  } else { plot.new(); title(expression(beta~"(t) unavailable")); box() }
  
  # (6) Det vs Stoch overlay (middle-right)
  if (!is.null(det_stoch) && exists("plot_det_vs_stoch", mode = "function")) {
    do.call(plot_det_vs_stoch, args = c(det_stoch, list()))
  } else { plot.new(); title("Det vs Stoch overlay (not provided)"); box() }
  
  # (7) Peak I
  if (has_det) {
    peak_I <- apply(sim$I, 2, max, na.rm = TRUE)
    bp <- barplot(peak_I, col = "gray80", border = "gray40", ylim = .safe_range(peak_I),
                  main = "Peak I by population", ylab = "peak I (proportion)")
    text(bp, peak_I, labels = sprintf("%.3f", peak_I), pos = 3, cex = 0.9)
  } else {
    I_arr <- sim$proportions[, , , "I"]; n_s <- dim(I_arr)[2]
    peak_mat <- apply(I_arr, c(2, 3), max, na.rm = TRUE)
    mean_p <- colMeans(peak_mat, na.rm = TRUE)
    qL_p <- apply(peak_mat, 2, .qfun, p = probs[1]); qU_p <- apply(peak_mat, 2, .qfun, p = probs[2])
    bp <- barplot(mean_p, col = "gray80", border = "gray40",
                  main = sprintf("Peak I by population (n=%d sims)", n_s),
                  ylab = sprintf("mean peak I (%.0f–%.0f%% band)", probs[1]*100, probs[2]*100),
                  ylim = .safe_range(0, qU_p))
    arrows(bp, qL_p, bp, qU_p, angle = 90, code = 3, length = 0.05, lwd = 2)
    text(bp, mean_p, labels = sprintf("%.3f", mean_p), pos = 3, cex = 0.9)
  }
  
  # (8) Stochasticity magnitude
  if (has_sto) {
    width_I <- colMeans(I_qU - I_qL, na.rm = TRUE)
    barplot(width_I, col = grDevices::adjustcolor("gray40", 0.6), border = "gray20",
            main = "Stochasticity (avg ribbon width of I)", ylab = "avg (qU - qL)",
            ylim = .safe_range(0, width_I))
  } else { plot.new(); title("Stochasticity (det only: none)"); box() }
  
  # (9) Parameters
  plot.new(); box(); title("Parameters", line = 0.5)
  usr <- par("usr"); xmid <- (usr[1] + usr[2]) / 2; y <- 0.92
  add_line <- function(lbl, val) { text(xmid, y, paste0(lbl, ": ", val), cex = 0.9); y <<- y - 0.095 }
  p <- sim$params
  add_line("n_times", length(time)); add_line("P", P)
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
