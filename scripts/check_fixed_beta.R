# --- Fixed-β sanity check: deterministic vs mean of stochastic ---

if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
R.utils::sourceDirectory("R/", modifiedOnly = FALSE)

# Source Christy's simulators (repo root)
source("simulate_sirs_det.R")
source("simulate_sirs_stoch.R")  # defines simulate_sirs(...)
# (β is constant here, so no need make_beta.R)

# ---- Parameters (keep simple & fast) ----
n_times <- 200
pop     <- 100000
I_init  <- 10
beta0   <- 0.18       # CONSTANT beta for sanity check
gamma   <- 1/7        # ~7-day infectious period
omega   <- 1/30       # ~30-day immunity duration
n_sims  <- 200        # number of stochastic runs
seed    <- 42

# ---- Run deterministic with constant β ----
det <- simulate_sirs_det(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta0, gamma = gamma, omega = omega
)

# ---- Run stochastic with same parameters & constant β ----
st <- simulate_sirs(
  n_times = n_times, pop = pop, I_init = I_init,
  beta = beta0, gamma = gamma, omega = omega,
  epsilon = 0, alpha = NULL,
  n_sims = n_sims, stochastic = TRUE, seed = seed
)

# ---- Helpers: mean & quantile bands across sims ----
col_mean <- function(mat_time_by_sims) rowMeans(mat_time_by_sims, na.rm = TRUE)
band <- function(mat_time_by_sims, probs = c(0.1, 0.9)) {
  qs <- apply(mat_time_by_sims, 1, quantile, probs = probs, na.rm = TRUE)
  list(low = qs[1, ], high = qs[2, ])
}
rmse <- function(a, b) sqrt(mean((a - b)^2))

# Extract I(t): proportions[time, sim, "I"] -> make [time x sims]
I_mat <- st$proportions[, , "I", drop = FALSE]
I_mat <- I_mat[, , 1]  # now [time x sims]

I_mean <- col_mean(I_mat)
I_band <- band(I_mat, probs = c(0.1, 0.9))

# Extract incidence(t): cases [time x sims]
Inc_mat  <- st$cases
Inc_mean <- col_mean(Inc_mat)
Inc_band <- band(Inc_mat, probs = c(0.1, 0.9))

# ---- Print sanity metrics ----
cat("RMSE (I proportion):     ", rmse(det$I, I_mean), "\n")
cat("RMSE (incidence counts): ", rmse(det$incidence, Inc_mean), "\n")

# ---- Plot: deterministic line over stochastic ribbon + mean ----
ribbon <- function(x, low, high, col = "gray", alpha = 0.25) {
  polygon(c(x, rev(x)), c(low, rev(high)), border = NA,
          col = adjustcolor(col, alpha.f = alpha))
}

op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# (1) I(t) proportions
ylim_I <- range(0, 1, det$I, I_band$high, na.rm = TRUE)
plot(det$time, det$I, type = "l", lwd = 2, col = "#0072B2",
     ylim = ylim_I, xlab = "day", ylab = "I (proportion)",
     main = "I(t): det vs stochastic")
ribbon(det$time, I_band$low, I_band$high)
lines(det$time, I_mean, lwd = 2, col = "black")   # stochastic mean
lines(det$time, det$I,   lwd = 2, col = "#0072B2")# deterministic
legend("topright", c("stoch 10-90%", "stoch mean", "deterministic"),
       lty = c(1,1,1), lwd = c(10,2,2),
       col = c(adjustcolor("gray",0.25),"black","#0072B2"),
       bty = "n", seg.len = 2)

# (2) Incidence counts
ylim_inc <- range(0, det$incidence, Inc_band$high, na.rm = TRUE)
plot(det$time, det$incidence, type = "l", lwd = 2, col = "#CC79A7",
     ylim = ylim_inc, xlab = "day", ylab = "new infections (count)",
     main = "Incidence: det vs stochastic")
ribbon(det$time, Inc_band$low, Inc_band$high)
lines(det$time, Inc_mean, lwd = 2, col = "black")   # stochastic mean
lines(det$time, det$incidence, lwd = 2, col = "#CC79A7")  # deterministic

# ---- Done ----
cat("\nFixed-β sanity check completed ✅\n")
