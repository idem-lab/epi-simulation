# Project: Kids Research Institute — SIRS modelling
# Script: scripts/demo_plot_det_vs_stoch.R
# Purpose: Visualise deterministic line over stochastic ribbon (I and incidence)
# Inputs: functions in R/ (auto-sourced) + Christy’s root simulators
# Outputs: plots/demo_det_vs_stoch.pdf

# Load all functions from R/
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
R.utils::sourceDirectory("R/", modifiedOnly = FALSE)

# Christy's simulators (still in repo root)
source("simulate_sirs_det.R")
source("simulate_sirs_stoch.R")

outdir <- "plots"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
pdf(file.path(outdir, "demo_det_vs_stoch.pdf"), width = 12, height = 5)
on.exit(dev.off(), add = TRUE)

# Params
n_times <- 200; pop <- 100000; I_init <- 10
beta0 <- 0.18; gamma <- 1/7; omega <- 1/30
n_sims <- 200; seed <- 42

# Deterministic (constant b