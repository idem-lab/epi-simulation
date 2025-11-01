# Project: Kids Research Institute — SIRS modelling
# Script: scripts/demo_plot_det_vs_stoch.R
# Output: plots/demo_det_vs_stoch.pdf

# ---- deps & paths ----
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
if (!requireNamespace("here", quietly = TRUE))     install.packages("here")
if (!requireNamespace("grDevices", quietly = TRUE)) stop("grDevices missing?")

library(R.utils); library(here)

# Source helpers from R/
R.utils::sourceDirectory(here("R"), modifiedOnly = FALSE)

# If these simulators are separate files, source them explicitly:
if (file.exists(here("R", "simulate_sirs_det.R")))   source(here("R", "simulate_sirs_det.R"))
if (file.exists(here("R", "simulate_sirs_stoch.R"))) source(here("R", "simulate_sirs_stoch.R"))

# ---- sanity checks ----
stopifnot(exists("simulate_sirs_det"))
stopifnot(exists("simulate_sirs"))         # or simulate_sirs_stoch — adjust name if needed
stopifnot(exists("plot_det_vs_stoch"))

outdir <- here("plots")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile <- file.path(outdir, "demo_det_vs_stoch.pdf")

# ---- params ----
n_times <- 200; pop <- 100000; I_init <- 10
beta0 <- 0.18; gamma <- 1/7; omega <- 1/30
n_sims <- 200; seed <- 42

# ---- run & plot ----
make_pdf <- function() {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outfile <- file.path(outdir, "demo_det_vs_stoch.pdf")
  
  # Prefer cairo_pdf; fall back to pdf()
  dev <- tryCatch(
    { grDevices::cairo_pdf(filename = outfile, width = 12, height = 5) },
    error = function(e) {
      message("cairo_pdf failed: ", e$message, " — falling back to pdf().")
      grDevices::pdf(file = outfile, width = 12, height = 5)
    }
  )
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    if (file.exists(outfile)) {
      cat("Wrote:", outfile, "size =", file.info(outfile)$size, "bytes\n")
    } else {
      cat("No output file found at:", outfile, "\n")
    }
  }, add = TRUE)
  
  # Simulate
  det <- simulate_sirs_det(n_times, pop, I_init, beta0, gamma, omega)
  st  <- simulate_sirs(n_times, pop, I_init, beta0, gamma, omega,
                       epsilon = 0, alpha = NULL, n_sims = n_sims,
                       stochastic = TRUE, seed = seed)
  
  # Plot
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(1, 2))
  plot_det_vs_stoch(det, st, state = "I",         probs = c(0.1, 0.9))
  plot_det_vs_stoch(det, st, state = "incidence", probs = c(0.1, 0.9))
}

make_pdf()
