# ================================================
# demo_run.R — SIRS demo (deterministic, stochastic, multi-pop)
# Sections:
#   A) Deterministic (constant beta)
#   B) Deterministic (seasonal beta)
#   C) Stochastic (seasonal beta) + per-sim summary
#   D) Sanity check (removed , seperate file)
#   E) Multi-pop deterministic
#   F) Multi-pop stochastic (mean + ribbons)
#   G) Dashboards (6-panel and 9-panel)
# ================================================

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

# Stochastic (C)
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

# ---- 2) RUN DEMO ----

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
  beta_mat  = mp_beta,     # scalar OK; recycled internally
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

# ---- DASHBOARDS ----
cat("\n[G] Dashboards\n")
# Build constant-β comparison for the overlay panel
det_const <- simulate_sirs_det(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta, gamma = gamma, omega = omega
)
st_const <- simulate_sirs(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta, gamma = gamma, omega = omega,
  epsilon = 0, alpha = NULL, n_sims = 200, stochastic = TRUE, seed = 123
)

# 6-panel dashboard
plot_dashboard(
  sim_E,
  probs = c(0.05, 0.95),
  per_million = TRUE,
  main = "Multi-pop deterministic",
  det_stoch = list(
    det     = det_const,
    stoch   = st_const,
    state   = "incidence",      # or "I"
    probs   = c(0.10, 0.90),
    det_col = "#2C7FB8",
    main    = "Det vs Stoch (constant β)"
  )
)

# 9-panel dashboard
plot_dashboard_v2(
  sim_F,
  probs = mp_ribbon,
  per_million = TRUE,
  main = "Multi-pop stochastic",
  det_stoch = list(
    det     = det_const,
    stoch   = st_const,
    state   = "incidence",      # overlay content; use "I" for prevalence
    probs   = c(0.10, 0.90),
    det_col = "#2C7FB8",
    main    = "Det vs Stoch (constant β)"
  )
)

# ==============================================================
# 3) INTERACTIVE PLOT PICKER — call run_plot_picker() from Console
# ==============================================================

plot_actions <- list(
  "A) Det. constant β — SIR both_side" = function() {
    plot_sir_diag(sim_A, which = "both_side")
  },
  "B) Det. seasonal β — overlay" = function() {
    plot_sir_diag(sim_B, which = "overlay")
  },
  "C) Stoch. seasonal β — first trajectory" = function() {
    plot_sir_diag(list(
      time = stoch_C$time,
      S = stoch_C$proportions[, 1, "S"],
      I = stoch_C$proportions[, 1, "I"],
      R = stoch_C$proportions[, 1, "R"],
      incidence = stoch_C$cases[, 1]
    ), which = "both_side")
  },
  "E1) Multi-pop det. — I(t) per pop" = function() {
    plot_multi_I(sim_E)
  },
  "E2) Multi-pop det. — incidence per million" = function() {
    plot_multi_incidence(sim_E, per_million = TRUE)
  },
  "F1) Multi-pop stoch. — ribbons for I(t)" = function() {
    plot_multi_stoch_I(sim_F, probs = mp_ribbon)
  },
  "F2) Multi-pop stoch. — ribbons for incidence per million" = function() {
    plot_multi_stoch_incidence(sim_F, per_million = TRUE, probs = mp_ribbon)
  },
  "G1) 6-panel dashboard" = function() {
    plot_dashboard(
      sim_E,
      probs = c(0.05, 0.95),
      per_million = TRUE,
      main = "Multi-pop deterministic",
      det_stoch = list(
        det     = det_const,
        stoch   = st_const,
        state   = "incidence",
        probs   = c(0.10, 0.90),
        det_col = "#2C7FB8",
        main    = "Det vs Stoch (constant β)"
      )
    )
  },
  "G2) 9-panel dashboard" = function() {
    plot_dashboard_v2(
      sim_F,
      probs = mp_ribbon,
      per_million = TRUE,
      main = "Multi-pop stochastic",
      det_stoch = list(
        det     = det_const,
        stoch   = st_const,
        state   = "incidence",
        probs   = c(0.10, 0.90),
        det_col = "#2C7FB8",
        main    = "Det vs Stoch (constant β)"
      )
    )
  }
)

run_plot_picker <- function() {
  if (!interactive()) stop("Interactive console required.")
  
  # 1) Try one-line multi-select first
  cat("\nSelect one or more plots to render:\n\n")
  for (i in seq_along(plot_actions)) cat(sprintf("%2d: %s\n", i, names(plot_actions)[i]))
  cat("\n")
  flush.console()
  ans <- tryCatch(readline("Enter numbers separated by space/comma (or press Enter to skip): "),
                  error = function(e) "")
  
  tokens <- unlist(strsplit(ans, "[,\\s]+"))
  nums <- unique(suppressWarnings(as.integer(tokens)))
  nums <- nums[!is.na(nums) & nums >= 1 & nums <= length(plot_actions)]
  
  if (length(nums) > 0) {
    cat("\nRendering:\n",
        paste(sprintf("%2d: %s", nums, names(plot_actions)[nums]), collapse = "\n"),
        "\n\n", sep = "")
    invisible(lapply(plot_actions[nums], function(f) f()))
    cat("\nAll done.\n")
    return(invisible(NULL))
  }
  
  # 2) Fallback: step-by-step menu loop (rock-solid)
  cat("\nNo line input detected. Switching to step-by-step selector.\n")
  picked <- integer(0)
  repeat {
    choice <- utils::menu(c(names(plot_actions), "Done"),
                          title = "Pick a plot (choose one at a time)")
    if (choice <= 0 || choice > length(plot_actions)) break
    picked <- c(picked, choice)
    plot_actions[[choice]]()  # draw immediately
  }
  if (!length(picked)) {
    cat("\nNo selections made. Nothing to render.\n")
  } else {
    cat("\nRendered:\n",
        paste(sprintf("%2d: %s", picked, names(plot_actions)[picked]), collapse = "\n"),
        "\n\n", sep = "")
  }
  cat("All done.\n")
}