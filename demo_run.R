source("simulate_sirs_det.R")
source("plot_sir_diag.R")
source("make_beta.R")


# A) deterministic, constant beta
sim1 <- simulate_sirs_det(
  n_times = 365, pop = 50000, I_init = 10,
  beta = 0.16, gamma = 1/20, omega = 1/7
)
plot_sir_diag(sim1, "both")

# B) deterministic, SEASONAL beta
beta2 <- make_beta(365, mode="seasonal", base=0.18, amplitude=0.25, phase=30)
sim2 <- simulate_sirs_det(
  n_times = 365, pop = 100000, I_init = 10,
  beta = beta2, gamma = 1/30, omega = 1/7
)
plot_sir_diag(sim2, "overlay")

source("simulate_sirs_stoch.R")
source("summarize_sim.R")

# C) stochastic: many runs, seasonal beta
beta3 <- make_beta(365, mode="seasonal", base=0.18, amplitude=0.25, phase=30)
st <- simulate_sirs(
  n_times=365, pop=100000, I_init=10,
  beta=beta3, gamma=1/7, omega=1/30,
  epsilon=1e-4, alpha=0.3,
  n_sims=100, stochastic=TRUE, seed=42
)

# Plot the first trajectory quickly using your existing plot helper
plot_sir_diag(list(
  time = st$time,
  S = st$proportions[,1,"S"],
  I = st$proportions[,1,"I"],
  R = st$proportions[,1,"R"],
  incidence = st$cases[,1]
), which="both")

# Summaries across sims
summary_tbl <- summarize_sim(st)
print(head(summary_tbl, 5))


