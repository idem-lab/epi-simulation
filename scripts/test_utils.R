# Project: Kids Research Institute — SIRS modelling
# Script: scripts/test_utils.R
# Purpose: Exercise Varun’s utilities across deterministic/stochastic/multi-pop
# Inputs: functions in R/ (auto-sourced) + Christy’s root simulators
# Outputs: console checks (optional plots to plots/ if added)

if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
R.utils::sourceDirectory("R/", modifiedOnly = FALSE)

# --- Source Christy's simulators (root) ---
source("simulate_sirs_det.R")
source("simulate_sirs_stoch.R")
source("simulate_sirs_multi_det.R")
source("make_beta.R")
source("make_contact.R")

# --- Source Varun's utilities (in R/) ---
source("R/reff_from_sim.R")
source("R/adjust_beta.R")
source("R/to_tidy.R")
source("R/check_contact.R")
source("R/cum_metrics.R")

cat("=== 1) Deterministic single-pop test ===\n")
# Build a seasonal beta and run sim
beta_det <- make_beta(365, mode = "seasonal", base = 0.18, amplitude = 0.25, phase = 30)
sim_det <- simulate_sirs_det(
  n_times = 365, pop = 100000, I_init = 10,
  beta = beta_det, gamma = 1/7, omega = 1/30
)

# Varun utils on deterministic output
reff_det <- reff_from_sim(sim_det)
df_det   <- to_tidy(sim_det)
ci_det   <- cumulative_incidence(sim_det)
ar_det   <- attack_rate(sim_det)

# Quick checks
stopifnot(nrow(reff_det) == length(sim_det$time))
stopifnot(all(c("time","group","sim","state","value") %in% names(df_det)))
stopifnot(is.numeric(ar_det), ar_det >= 0, ar_det <= 1)

# Simple scenario: reduce beta by 30% for days 60–88
beta_det_2 <- adjust_beta(beta_det, data.frame(start = 60, end = 88), multiplier = 0.7)
sim_det_2 <- simulate_sirs_det(365, 100000, 10, beta_det_2, 1/7, 1/30)

cat("Deterministic peak incidence (baseline vs scenario):",
    max(sim_det$incidence), "->", max(sim_det_2$incidence), "\n")

cat("=== 2) Stochastic multi-run test ===\n")
# Use 180 days to keep it quick
beta_st <- make_beta(180, mode = "seasonal", base = 0.18, amplitude = 0.25, phase = 30)
st <- simulate_sirs(
  n_times = 180, pop = 100000, I_init = 10,
  beta = beta_st, gamma = 1/7, omega = 1/30,
  epsilon = 0, alpha = NULL, n_sims = 5, stochastic = TRUE, seed = 42
)

reff_st <- reff_from_sim(st)
df_st   <- to_tidy(st)
ci_st   <- cumulative_incidence(st)
ar_st   <- attack_rate(st)

# Checks
stopifnot(all(c("time","sim","Reff") %in% names(reff_st)))
stopifnot(length(ar_st) == 5)

cat("Stochastic: attack rate per simulation (last-day R):\n")
print(round(ar_st, 3))

cat("=== 3) Multi-population deterministic test ===\n")
# 3 groups with mixing
C <- make_contact(P = 3, within = 0.9)
check_contact(C)                         # validate
C_row <- normalize_contact(C, "row")     # normalized (optional)

beta_mat <- matrix(0.20, nrow = 200, ncol = 3)
# Scenario: reduce beta in all groups on days 50–80
beta_mat2 <- adjust_beta(beta_mat, data.frame(start = 50, end = 80), multiplier = 0.75)

sim_multi <- simulate_sirs_multi(
  n_times = 200,
  pop_vec = c(100000, 50000, 20000),
  I_init  = c(10, 5, 0),
  beta_mat = beta_mat2,
  gamma = 1/7, omega = 1/30,
  C = C_row
)

reff_multi <- reff_from_sim(sim_multi)
df_multi   <- to_tidy(sim_multi)
ci_multi   <- cumulative_incidence(sim_multi)
ar_multi   <- attack_rate(sim_multi)

# Checks
stopifnot(all(c("time","group","Reff") %in% names(reff_multi)))
stopifnot(length(ar_multi) == 3)

cat("Multi-pop: final attack rate by group:\n")
print(round(ar_multi, 3))

cat("\nAll Varun utility tests completed successfully ✅\n")
