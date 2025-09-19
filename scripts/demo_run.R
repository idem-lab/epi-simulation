#source("simulate_sirs_det.R")
#source("plot_sir_diag.R")
#source("make_beta.R")


# A) deterministic, constant beta
#sim1 <- simulate_sirs_det(
#  n_times = 365, pop = 50000, I_init = 10,
#  beta = 0.16, gamma = 1/20, omega = 1/7
#)
#plot_sir_diag(sim1, "both")

# B) deterministic, SEASONAL beta
#beta2 <- make_beta(365, mode="seasonal", base=0.18, amplitude=0.25, phase=30)
#sim2 <- simulate_sirs_det(
#  n_times = 365, pop = 100000, I_init = 10,
#  beta = beta2, gamma = 1/30, omega = 1/7
#)
#plot_sir_diag(sim2, "overlay")

#source("simulate_sirs_stoch.R")
#source("summarize_sim.R")

# C) stochastic: many runs, seasonal beta
#beta3 <- make_beta(365, mode="seasonal", base=0.18, amplitude=0.25, phase=30)
#st <- simulate_sirs(
#  n_times=365, pop=100000, I_init=10,
#  beta=beta3, gamma=1/7, omega=1/30,
#  epsilon=1e-4, alpha=0.3,
#  n_sims=100, stochastic=TRUE, seed=42
#)

# Plot the first trajectory quickly using your existing plot helper
#plot_sir_diag(list(
#  time = st$time,
#  S = st$proportions[,1,"S"],
#  I = st$proportions[,1,"I"],
#  R = st$proportions[,1,"R"],
#  incidence = st$cases[,1]
#), which="both")

# Summaries across sims
#summary_tbl <- summarize_sim(st)
#print(head(summary_tbl, 5))



# ================================================
# Sections:
#   A) Deterministic (constant beta)
#   B) Deterministic (seasonal beta)
#   C) Stochastic (seasonal beta) + per-sim summary
#   D) Sanity check (removed , seperate file)
#   E) Multi-pop deterministic
#   F) Multi-pop stochastic (mean + ribbons)

# ---- 1. EDIT YOUR INPUTS HERE ----
# Common (single-pop)
n_times <- 365
pop     <- 100000
I_init  <- 10
beta    <- 0.16
gamma   <- 1/7
omega   <- 1/30

# Seasonal beta (for B, C)
season_base  <- 0.18
season_amp   <- 0.25     # 0 ≤ amplitude < 1
season_phase <- 30

# Stochastic (C, D)
n_sims   <- 200
epsilon  <- 1e-4
alpha    <- 0.3          
seed     <- 42

# Multi-pop (E, F)
mp_n_times <- 365
mp_pops    <- c(100000, 500, 2000)
mp_Iinit   <- c(10, 5, 0)
mp_beta    <- 0.20       # scalar; recycled to [time x pops]
mp_gamma   <- 1/7
mp_omega   <- 1/30
mp_contact <- matrix(c(1,0.2,0.1,
                       0.2,1,0.3,
                       0.1,0.3,1),
                     3, 3, byrow = TRUE)

# Multi-pop stochastic extra
mp_sims    <- 200
mp_epsilon <- 0
mp_alpha   <- NULL
mp_seed    <- 99
mp_ribbon  <- c(0.05, 0.95)   

# ---- 2) LOAD FUNCTIONS ----
source("simulate_sirs_det.R")
source("simulate_sirs_stoch.R")        
source("simulate_sirs_multi_det.R")
source("simulate_sirs_multi_stoch.R")  
source("make_beta.R")
source("plot_sir_diag.R")
source("plot_multi.R")                 
source("summarize_sim.R")
source("plot_dashboard.R")

# ---- 3) RUN DEMO ----

# A) Deterministic (constant beta)
cat("\n[A] Deterministic (constant beta)\n")
sim_A <- simulate_sirs_det(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta, gamma = gamma, omega = omega
)
plot_sir_diag(sim_A, which = "both_side")

# B) Deterministic (seasonal beta)
cat("\n[B] Deterministic (seasonal beta)\n")
beta_B <- make_beta(
  n_times, mode = "seasonal",
  base = season_base, amplitude = season_amp, phase = season_phase
)
sim_B <- simulate_sirs_det(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta_B, gamma = gamma, omega = omega
)
plot_sir_diag(sim_B, which = "overlay")

# C) Stochastic (seasonal beta) + per-sim summaries
cat("\n[C] Stochastic (seasonal beta)\n")
beta_C <- make_beta(
  n_times, mode = "seasonal",
  base = season_base, amplitude = season_amp, phase = season_phase
)
stoch_C <- simulate_sirs(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta_C, gamma = gamma, omega = omega,
  epsilon = epsilon, alpha = alpha,
  n_sims = n_sims, stochastic = TRUE, seed = seed
)
# quick view: first trajectory
plot_sir_diag(list(
  time = stoch_C$time,
  S = stoch_C$proportions[, 1, "S"],
  I = stoch_C$proportions[, 1, "I"],
  R = stoch_C$proportions[, 1, "R"],
  incidence = stoch_C$cases[, 1]
), which = "both_side")
# headline summaries (first few sims)
print(utils::head(summarize_sim(stoch_C), 5))

# E) Multi-pop deterministic
cat("\n[E] Multi-pop deterministic\n")
beta_mat <- matrix(mp_beta, nrow = mp_n_times, ncol = length(mp_pops))
sim_E <- simulate_sirs_multi(
  n_times = mp_n_times, pop_vec = mp_pops, I_init = mp_Iinit,
  beta_mat = beta_mat, gamma = mp_gamma, omega = mp_omega, C = mp_contact
)
plot_multi_I(sim_E)
plot_multi_incidence(sim_E, per_million = TRUE)

# F) Multi-pop stochastic (mean + ribbons)
cat("\n[F] Multi-pop stochastic\n")
sim_F <- simulate_sirs_multi_stoch(
  n_times   = mp_n_times,
  pop_vec   = mp_pops,
  I_init    = mp_Iinit,
  beta_mat  = mp_beta,     
  gamma     = mp_gamma,
  omega     = mp_omega,
  C         = mp_contact,
  epsilon   = mp_epsilon,
  alpha     = mp_alpha,
  n_sims    = mp_sims,
  stochastic = TRUE,
  seed      = mp_seed
)
# ribbons for I(t) per population
plot_multi_stoch_I(sim_F, probs = mp_ribbon)
# ribbons for daily incidence per population (per million)
plot_multi_stoch_incidence(sim_F, per_million = TRUE, probs = mp_ribbon)


# ----DASHBOARDS ----
# Deterministic dashboard
plot_dashboard(sim_E, probs = c(0.05, 0.95), per_million = TRUE,
               main = "Multi-pop deterministic")

# Stochastic dashboard
plot_dashboard(sim_F, probs = mp_ribbon, per_million = TRUE,
               main = "Multi-pop stochastic")


plot_dashboard_v2(sim_E, probs = c(0.05, 0.95), per_million = TRUE,
                  main = "Multi-pop deterministic")


plot_dashboard_v2(sim_F, probs = mp_ribbon, per_million = TRUE,
                  main = "Multi-pop stochastic")