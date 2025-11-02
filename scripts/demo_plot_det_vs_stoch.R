# =========================================
# Det vs Stoch SIRS — tight visual match
# =========================================

# ---- Controls you might tweak ----
dt           <- 0.1
n_sims       <- 2000
ribbon_probs <- c(0.35, 0.65)   # narrower band
center       <- "median"        # or "mean"
beta_nudge   <- 1.0             # e.g., 0.985 for tiny presentation nudge

# ---- Model params ----
n_times <- 365
pop     <- 100000
I_init  <- 10
beta    <- 0.18 * beta_nudge
gamma   <- 1/7
omega   <- 1/30

# ---- Discrete-time, rate->prob mapping helpers ----
simulate_sirs_det_discrete <- function(n_times, pop, I_init, beta, gamma, omega, dt = 0.25) {
  n_steps <- ceiling(n_times / dt)
  S <- numeric(n_steps + 1); I <- S; R <- S
  inc <- numeric(n_steps)
  
  S[1] <- pop - I_init; I[1] <- I_init; R[1] <- 0
  
  for (k in 1:n_steps) {
    p_inf <- 1 - exp(-(beta * I[k] / pop) * dt)
    p_rec <- 1 - exp(-gamma * dt)
    p_wan <- 1 - exp(-omega * dt)
    
    new_inf <- S[k] * p_inf
    new_rec <- I[k] * p_rec
    new_wan <- R[k] * p_wan
    
    S[k+1] <- S[k] - new_inf + new_wan
    I[k+1] <- I[k] + new_inf - new_rec
    R[k+1] <- R[k] + new_rec - new_wan
    
    inc[k] <- new_inf
  }
  
  spd <- round(1/dt)
  days <- 1:floor(n_steps/spd)
  daily_inc <- tapply(inc[1:(tail(days,1)*spd)], rep(days, each = spd), sum)
  
  list(tim