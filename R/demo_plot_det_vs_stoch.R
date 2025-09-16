# Demo for plot_det_vs_stoch

# Christy's simulators
source("simulate_sirs_det.R")
source("simulate_sirs_stoch.R")
source("make_beta.R")

# Varun's plot
source("R/plot_det_vs_stoch.R")

# Params
n_times <- 200; pop <- 100000; I_init <- 10
beta0 <- 0.18; gamma <- 1/7; omega <- 1/30
n_sims <- 200; seed <- 42

# Deterministic (constant beta for clean comparison)
det <- simulate_sirs_det(n_times, pop, I_init, beta0, gamma, omega)

# Stochastic with same params
st <- simulate_sirs(n_times, pop, I_init, beta0, gamma, omega,
                    epsilon = 0, alpha = NULL, n_sims = n_sims,
                    stochastic = TRUE, seed = seed)

# Two panels for I and incidence
op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
par(mfrow = c(1,2))
plot_det_vs_stoch(det, st, state = "I",         probs = c(0.1, 0.9))
plot_det_vs_stoch(det, st, state = "incidence", probs = c(0.1, 0.9))
