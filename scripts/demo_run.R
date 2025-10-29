# Sections:
#(A) Deterministic (constant beta)
#(B) Deterministic (seasonal beta)
#(C) Stochastic (constant beta)
#(D) Stochastic (seasonal beta)
#(E) Multi-pop deterministic
#(F) Multi-pop stochastic
#(G) Summary table of simulation results
#(H) Dashboards with self arrangement

# ---- EDIT YOUR INPUTS HERE -------------------
# Common components
n_times <- 60      # Total days to simulate
pop     <- 100000  # Population size
I_init  <- 10      # Initial number of infectious individuals
beta    <- 0.75    # Transmission rate  
gamma   <- 1/7     # Recovery rate
omega   <- 1/30    # Rate of loss of immunity(R -> S)


# Seasonal beta (for B, D)
season_base  <- 0.70   # Base transmission rate
season_amp   <- 0.25   # Amplitude of seasonal variation (0 ≤ amplitude < 1)
season_phase <- 30     # Phase shift in days

# Stochastic (C,D,F)
n_sims   <- 200   # Number of independent simulation runs (columns).
epsilon  <- 1e-4  # External infection pressure.
alpha    <- 0.3   # Reporting rate
seed     <- 42    # Random seed for reproducibility

# Multi-pop (E, F)
mp_n_times <- 180                     # Total days to simulate
mp_pops    <- c(100000, 500, 200000)  # Population sizes
mp_Iinit   <- c(100, 5, 1000)         # Initial infectious individuals
mp_beta    <- 0.20                    # Transmission rate
mp_gamma   <- 1/7                     # Recovery rate
mp_omega   <- 1/30                    # Rate of loss of immunity(R -> S)

# Multi-pop stochastic (F)
mp_sims    <- 200      # Number of stochastic simulations
mp_epsilon <- 0        # Noise parameter
mp_alpha   <- NULL     # Scaling factor for stochasticity
mp_seed    <- 99       # Random seed for reproducibility

# Customise uncertainty ribbons(Default to 95%)
ribbon_probs <- c(0.025, 0.975)

# Customise colors for multiple populations
my_pop_cols <- c("#E41A1C", "#377EB8", "#4DAF4A")

# ---- 2) LOAD FUNCTIONS --------------------------
R.utils::sourceDirectory("R/")

# ---- 3) RUN DEMO --------------------------------

# A) Deterministic (constant beta)
cat("\n[A] Deterministic (constant beta)\n")
sim_A <- simulate_sirs_det(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta, gamma = gamma, omega = omega
)
plot_sirs(sim_A, which = "overlay")
plot_sirs(sim_A, which = "sir")
plot_sirs(sim_A, which = "incidence")
plot_sirs(sim_A, which = "both_side")

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
plot_sirs(sim_B, which = "overlay")
plot_sirs(sim_B, which = "sir")
plot_sirs(sim_B, which = "incidence")
plot_sirs(sim_B, which = "both_side")

# C) Stochastic (constant beta) — single pop
cat("\n[C] Stochastic (constant beta)\n")
stoch_C <- simulate_sirs_stoch(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta, gamma = gamma, omega = omega,
  epsilon = epsilon, alpha = alpha,
  n_sims = n_sims, stochastic = TRUE, seed = seed + 1
)
stoch_C$params$ribbon_probs <- ribbon_probs

plot_stoch(stoch_C, which = "SIR")
plot_stoch(stoch_C, which = "incidence", per_million = TRUE)
plot_stoch(stoch_C, which = "overlay")
plot_stoch(stoch_C, which = "both")

# D) Stochastic (seasonal beta) — single pop
cat("\n[D] Stochastic (seasonal beta)\n")
beta_CS <- make_beta(
  n_times, mode = "seasonal",
  base = season_base, amplitude = season_amp, phase = season_phase
)
stoch_CS <- simulate_sirs_stoch(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta_CS, gamma = gamma, omega = omega,
  epsilon = epsilon, alpha = alpha,
  n_sims = n_sims, stochastic = TRUE, seed = seed
)
stoch_CS$params$ribbon_probs <- ribbon_probs

plot_stoch(stoch_CS, which = "incidence")
plot_stoch(stoch_CS, which = "SIR")
plot_stoch(stoch_CS, which = "overlay")
plot_stoch(stoch_CS, which = "both")

# E) Multi-pop deterministic
cat("\n[E] Multi-pop deterministic\n")
beta_mat <- matrix(mp_beta, nrow = mp_n_times, ncol = length(mp_pops))
sim_E <- simulate_sirs_multi(
  n_times = mp_n_times, pop_vec = mp_pops, I_init = mp_Iinit,
  beta_mat = beta_mat, gamma = mp_gamma, omega = mp_omega
)

plot_multi(sim_E, which = "S", group_style = "combined")
plot_multi(sim_E, which = "I", group_style = "facet")

# F) Multi-pop stochastic
cat("\n[F] Multi-pop stochastic\n")
sim_F <- simulate_sirs_multi_stoch(
  n_times    = mp_n_times,
  pop_vec    = mp_pops,
  I_init     = mp_Iinit,
  beta_mat   = mp_beta,     
  gamma      = mp_gamma,
  omega      = mp_omega,
  epsilon    = mp_epsilon,
  alpha      = mp_alpha,
  n_sims     = mp_sims,
  stochastic = TRUE,
  seed       = mp_seed
)
sim_F$params$ribbon_probs <- ribbon_probs

plot_multi(sim_F, which = "R", group_style = "combined")
plot_multi(sim_F, which = "incidence", group_style = "facet", per_million = TRUE)


# G) Dashboards — multi-pop
cat("\n[G] Dashboards — multi-pop \n")
plots_E <- plot_dashboard(
  sim           = sim_E,
  probs         = ribbon_probs, 
  per_million   = TRUE,
  group_style   = "combined",
  pop_cols      = my_pop_cols,
  show_bands    = FALSE                 
)
# Keep only implemented panels
panels_E <- c("S","I","R","incidence","beta","params")
p_E <- arrange_dashboard(plots_E[panels_E], layout = c(3,2), collect_legend = TRUE)
print(p_E)

panels_E2 <- c("S","I","R","incidence")
p_E2 <- arrange_dashboard(plots_E[panels_E2], layout = c(2,2), collect_legend = TRUE)
print(p_E2)

plots_F <- plot_dashboard(
  sim           = sim_F,
  probs         = ribbon_probs,
  per_million   = TRUE,
  group_style   = "facet",
  pop_cols      = my_pop_cols,
  show_bands    = TRUE                    # stochastic ribbons visible
)
# Keep only implemented panels
panels_F <- c("S","I","R","incidence","beta","params")
p_F <- arrange_dashboard(plots_F[panels_F], layout = c(3,2), collect_legend = TRUE)
print(p_F)

panels_F2 <- c("S","I","R","incidence")
p_F2 <- arrange_dashboard(plots_F[panels_F2], layout = c(2,2), collect_legend = TRUE)
print(p_F2)

# H) Single-pop “basic SIR” dashboard example
cat("\n[H] Single-pop dashboard — basic SIR panel (lines only)\n")
plots_single <- plot_dashboard(
  sim           = sim_B,
  probs         = ribbon_probs,
  per_million   = FALSE,
  group_style   = "facet",
  show_bands    = FALSE                   # deterministic: lines only
)
panels_single <- c("sir_basic","incidence","beta","params")
p_single <- arrange_dashboard(plots_single[panels_single], layout = c(2,2), collect_legend = TRUE)
print(p_single)

cat("\nDone.\n")


stoch_G<- simulate_sirs_stoch(
  n_times = 180,pop = 100, I_init = 2,
  beta = 0.2, gamma = 1/30, omega = 1/14,
  epsilon = 0, alpha = 0.1,
  n_sims = n_sims, stochastic = TRUE, seed = seed + 1
)
stoch_G$params$ribbon_probs <- ribbon_probs
plot_stoch(stoch_G, which = "SIR")