source("simulate_sirs_multi_det.R")
source("make_contact.R")
source("plot_multi.R")

dim(sim$I)
sim <- simulate_sirs_multi(
  n_times = 365,
  pop_vec = c(1e5, 5e4, 2e4),
  I_init = c(10, 5, 0),
  beta_mat = matrix(0.2, 365, 3),
  gamma = 1/7,
  omega = 1/30,
  C = matrix(c(1,0.2,0.1,
               0.2,1,0.3,
               0.1,0.3,1), 3, 3, byrow=TRUE)
)
plot_multi_I(sim)


sim2 <- simulate_sirs_multi(
  n_times = 200,
  pop_vec = c(1e5, 5e4),
  I_init  = c(50, 0),                 # infection starts only in pop1
  beta_mat = matrix(0.18, 200, 2),
  gamma = 1/7,
  omega = 1/30,
  C = matrix(c(1, 0.05,   # Pop1 mostly mixes within itself
               0.05, 1),  # Pop2 also
             2, 2, byrow = TRUE)
)
plot_multi_I(sim2)

sim3 <- simulate_sirs_multi(
  n_times = 300,
  pop_vec = c(3e4, 7e4),              # children = 30k, adults = 70k
  I_init  = c(20, 5),                 # more cases among children
  beta_mat = matrix(0.25, 300, 2),
  gamma = 1/5,                        # faster recovery (5 days)
  omega = 1/90,                       # slower waning
  C = matrix(c(1.0, 0.4,
               0.4, 0.8), 2, 2, byrow = TRUE)
)
plot_multi_I(sim3)

