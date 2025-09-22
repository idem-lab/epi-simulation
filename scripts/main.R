# Project: Kids Research Institute — SIRS modelling
# Script: scripts/main.R
# Purpose: End-to-end run — simulate, sanity check, and quick plots
# Inputs: functions in R/ (auto-sourced)
# Outputs: plots/main_run.pdf + TRUE/FALSE sanity check

if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
R.utils::sourceDirectory("R/", modifiedOnly = FALSE)

# Params
set.seed(12)
n_times <- 365
pop     <- 100000
I_init  <- 10
omega   <- 1/30
gamma   <- 1/7

# Simple beta (lognormal) like your original
beta_vec <- rlnorm(n_times, meanlog = 1.6, sdlog = 0.03) / 10

# Run
sim <- simulate_sirs(
  n_times = n_times, pop = pop, I_init = I_init,
  omega = omega, gamma = gamma, beta_vec = beta_vec, seed = 12
)

# Sanity
print(sanity_check(sim))

outdir <- "plots"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
pdf(file.path(outdir, "main_run.pdf"), width = 8, height = 6)
on.exit(dev.off(), add = TRUE)

# Quick plots (base R)
plot(sim$t, sim$N_t, type = "l", cex = 0.5, pch = 20,
     ylab = "count", xlab = "time", main = "incident infections")

plot(sim$t, sim$S, type = "l", cex = 0.5, pch = 20, col = "#009E73",
     ylab = "proportion", xlab = "time")
points(sim$t, sim$I, type = "l", cex = 0.5, pch = 20, col = "#CC79A7")
points(sim$t, sim$R, type = "l", cex = 0.5, pch = 20, col = "#0072B2")
legend('topright', c("Susceptible","Infected","Recovered"),
       lty = 1, pch = 20, col = c("#009E73","#CC79A7","#0072B2"),
       title = "Proportion", cex = 1)
