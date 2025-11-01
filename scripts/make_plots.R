# scripts/make_plots.R — publication-ready PNGs only
# ------------------------------------------------------------
message(">> Sourcing functions from R/ ...")
R.utils::sourceDirectory("R/", modifiedOnly = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------- Visual helpers ----------
theme_epi <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25),
      panel.grid.major.y = element_line(linewidth = 0.25),
      plot.title        = element_text(face = "bold"),
      plot.subtitle     = element_text(margin = margin(b = 6)),
      plot.caption      = element_text(size = rel(0.85)),
      legend.position   = "bottom",
      legend.title      = element_blank()
    )
}

# colors for S, I, R and populations
col_SIR  <- c(S = "#1f78b4", I = "#e377c2", R = "#2ca02c")
pop_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF")

# Save helper (captures return-or-draw plots)
save_png <- function(expr, file, width=7.8, height=4.8, dpi=320, bg="white"){
  if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE)
  grDevices::png(file, width = width*dpi, height = height*dpi, res = dpi, bg = bg)
  on.exit(grDevices::dev.off(), add = TRUE)
  obj <- try(eval(expr), silent = TRUE)
  if (inherits(obj, "ggplot") || inherits(obj, "patchwork") ||
      inherits(obj, "grob")   || inherits(obj, "gTree") ||
      inherits(obj, "gList")  || inherits(obj, "gtable")) {
    try(print(obj), silent = TRUE)
  }
  message("   ✓ ", file)
}

annotate_peak <- function(x, y){
  i <- which.max(y)
  list(
    geom_vline(xintercept = x[i], linetype = 3, linewidth = 0.4),
    annotate("label", x = x[i], y = y[i], label = paste0("Peak: day ", x[i]),
             vjust = -0.6, size = 3, label.size = 0)
  )
}

# ---------- Deterministic single-pop ----------
n_times <- 60; pop <- 100000; I_init <- 10; beta <- 0.75; gamma <- 1/7; omega <- 1/30
simA <- simulate_sirs_det(n_times, pop, I_init, beta, gamma, omega)

# Build tidy SIR from simA (assumes proportions)
sirA <- tibble(
  day = simA$time %||% seq_len(n_times),
  S = simA$S %||% simA$proportions[,1,"S"],
  I = simA$I %||% simA$proportions[,1,"I"],
  R = simA$R %||% simA$proportions[,1,"R"]
) |>
  pivot_longer(-day, names_to = "state", values_to = "value")

# Incidence (counts) from ΔI + γ I_{t-1}
Ivec <- sirA |> filter(state == "I") |> arrange(day) |> pull(value)
Ilag <- dplyr::lag(Ivec, default = Ivec[1])
inc_counts <- (Ivec - Ilag) + gamma * Ilag
inc_counts <- pmax(inc_counts, 0) * pop
incA <- tibble(day = unique(sirA$day), inc = inc_counts)

# SIR only
save_png(quote(
  ggplot(sirA, aes(day, value, colour = state)) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = col_SIR) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = "S, I, R", x = "Day", y = "Proportion") +
    theme_epi()
), "plots/det_sir.png")

# Incidence only
save_png(quote(
  ggplot(incA, aes(day, inc)) +
    geom_line(linewidth = 1, colour = "grey40") +
    annotate_peak(incA$day, incA$inc) +
    labs(title = "Daily incidence", x = "Day", y = "New infections (count)") +
    theme_epi()
), "plots/det_incidence.png")

# 2-panel composite (SIR + Incidence)
save_png(quote({
  p1 <- ggplot(sirA, aes(day, value, colour = state)) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = col_SIR) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = "S, I, R", x = "Day", y = "Proportion") +
    theme_epi()
  p2 <- ggplot(incA, aes(day, inc)) +
    geom_line(linewidth = 1, colour = "grey40") +
    annotate_peak(incA$day, incA$inc) +
    labs(title = "Daily incidence", x = "Day", y = "New infections (count)") +
    theme_epi()
  p1 + p2
}), "plots/det_both.png")

# ---------- Deterministic with seasonal beta ----------
season_base <- 0.70; season_amp <- 0.25; season_phase <- 30
beta_B <- make_beta(n_times, mode = "seasonal", base = season_base, amplitude = season_amp, phase = season_phase)
simB <- simulate_sirs_det(n_times, pop, I_init, beta_B, gamma, omega)
sirB <- tibble(
  day = simB$time %||% seq_len(n_times),
  S = simB$S %||% simB$proportions[,1,"S"],
  I = simB$I %||% simB$proportions[,1,"I"],
  R = simB$R %||% simB$proportions[,1,"R"]
) |> pivot_longer(-day, names_to = "state", values_to = "value")

save_png(quote(
  ggplot(sirB, aes(day, value, colour = state)) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = col_SIR) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Seasonal SIRS — S, I, R",
      subtitle = sprintf("β0=%.2f, amp=%.2f, phase=%d days", season_base, season_amp, season_phase),
      x = "Day", y = "Proportion"
    ) +
    theme_epi()
), "plots/seasonal_sir.png")

# ---------- Stochastic single-pop ----------
n_sims <- 200; epsilon <- 0; alpha <- 0.1; seed <- 42
stochC <- simulate_sirs(
  n_times = 80, pop = 1000, I_init = 5,
  beta = 0.2, gamma = 1/30, omega = 1/14,
  epsilon = epsilon, alpha = alpha,
  n_sims = n_sims, stochastic = TRUE, seed = seed
)
stochC$params$ribbon_probs <- c(0.025, 0.975)
tC <- stochC$time
probs <- stochC$params$ribbon_probs

# I(t) proportion ribbon
Imat <- stochC$proportions[, , "I", drop = TRUE]
I_med <- apply(Imat, 1, median, na.rm = TRUE)
I_lo  <- apply(Imat, 1, quantile, probs[1], na.rm = TRUE)
I_hi  <- apply(Imat, 1, quantile, probs[2], na.rm = TRUE)
dfI <- tibble(day = tC, med = I_med, lo = I_lo, hi = I_hi)

save_png(quote(
  ggplot(dfI, aes(day, med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, fill = "grey60") +
    geom_line(linewidth = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = "Stochastic — I(t) proportion", subtitle = "Median with 95% ribbon",
         x = "Day", y = "I proportion") +
    theme_epi()
), "plots/stoch_sir.png")

# Daily cases ribbon (reported)
Cmat <- stochC$cases
C_med <- apply(Cmat, 1, median, na.rm = TRUE)
C_lo  <- apply(Cmat, 1, quantile, probs[1], na.rm = TRUE)
C_hi  <- apply(Cmat, 1, quantile, probs[2], na.rm = TRUE)
dfC <- tibble(day = tC, med = C_med, lo = C_lo, hi = C_hi)

save_png(quote(
  ggplot(dfC, aes(day, med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, fill = "grey60") +
    geom_line(linewidth = 1) +
    labs(title = "Stochastic — daily incidence", subtitle = "Median reported cases with 95% ribbon",
         x = "Day", y = "Daily cases (reported)") +
    theme_epi()
), "plots/stoch_incidence.png")

if (exists("to_tidy", mode = "function")) {
  try(utils::write.csv(to_tidy(stochC), "plots/stochC_tidy_sample.csv", row.names = FALSE), silent = TRUE)
}

# ---------- Multi-pop deterministic ----------
mp_n_times <- 180
mp_pops    <- c(100000, 500, 200000)
mp_Iinit   <- c(100, 5, 1000)
mp_beta    <- 0.20
mp_gamma   <- 1/7
mp_omega   <- 1/30
beta_md <- matrix(mp_beta, nrow = mp_n_times, ncol = length(mp_pops))

simE <- simulate_sirs_multi(
  n_times = mp_n_times, pop_vec = mp_pops, I_init = mp_Iinit,
  beta_mat = beta_md, gamma = mp_gamma, omega = mp_omega
)

# Tidy S/I/R per group (force group to factor for discrete scales)
get_multi_states <- function(sim){
  if (exists("to_tidy", mode = "function")) {
    to_tidy(sim) |>
      filter(state %in% c("S","I","R")) |>
      mutate(group = as.factor(group))
  } else {
    t <- seq_len(nrow(sim$I)); groups <- seq_len(ncol(sim$I))
    bind_rows(lapply(groups, function(g){
      tibble(time = t, group = factor(g),
             S = sim$S[,g], I = sim$I[,g], R = sim$R[,g])
    })) |>
      pivot_longer(cols = c(S,I,R), names_to = "state", values_to = "value") |>
      mutate(group = as.factor(group))
  }
}
tdE <- get_multi_states(simE) |> mutate(group = as.factor(group))

# Helper for consistent manual palette keyed by factor levels
pal_for <- function(fct) {
  lv <- levels(fct)
  setNames(pop_cols[seq_len(length(lv))], lv)
}

# S (combined)
save_png(quote(
  tdE |>
    dplyr::filter(state == "S") |>
    ggplot(aes(time, value, colour = group)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = pal_for(tdE$group)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Multi-population — S (combined)", x = "Day", y = "S proportion", colour = "Group") +
    theme_epi()
), "plots/multi_S_combined.png")

# I (facet)
save_png(quote(
  tdE |>
    dplyr::filter(state == "I") |>
    ggplot(aes(time, value, colour = group)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = pal_for(tdE$group)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    facet_wrap(~group, scales = "free_y") +
    labs(title = "Multi-population — I (facet)", x = "Day", y = "I proportion", colour = "Group") +
    theme_epi()
), "plots/multi_I_facet.png")

# Derived incidence per million (ΔI + γ I_{t−1}), using factor-safe pop indexing
wideE <- tdE |>
  dplyr::filter(state == "I") |>
  dplyr::select(time, group, value) |>
  dplyr::group_by(group) |>
  dplyr::arrange(time, .by_group = TRUE) |>
  dplyr::mutate(Ilag = dplyr::lag(value, default = dplyr::first(value)),
                inc  = (value - Ilag) + mp_gamma * Ilag) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    pop = mp_pops[as.integer(group)],
    inc_per_million = inc / pop * 1e6
  )

save_png(quote(
  ggplot(wideE, aes(time, inc_per_million, colour = group)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = pal_for(wideE$group)) +
    labs(title = "Multi-population deterministic — incidence per million (derived)",
         x = "Day", y = "Incidence per million", colour = "Group") +
    theme_epi()
), "plots/multi_incidence_combined.png")

message(">> Done. PNGs saved under plots/.")
