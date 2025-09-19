simulate_sirs_multi_stoch <- function(
    n_times   = 365,
    pop_vec   = c(5e4, 5e4),
    I_init    = c(10, 10),
    beta_mat  = NULL,        # scalar | vector(n_times) | matrix[n_times x P]
    gamma     = 1/7,
    omega     = 1/30,
    C         = NULL,        # [P x P], rows=receivers, cols=sources
    epsilon   = 0,           # scalar or length-P
    alpha     = NULL,        # NULL or scalar/length-P (for cases only)
    n_sims    = 100,
    stochastic = TRUE,
    seed      = NULL
) {
  P <- length(pop_vec)
  stopifnot(n_times >= 2, P >= 1, length(I_init) == P)
  stopifnot(all(I_init >= 0 & I_init <= pop_vec))
  stopifnot(is.numeric(gamma) && gamma >= 0 && gamma <= 1)
  stopifnot(is.numeric(omega) && omega >= 0 && omega <= 1)
  if (!is.null(seed)) set.seed(seed)
  
  # Contact matrix
  if (is.null(C)) C <- matrix(1, nrow = P, ncol = P)
  stopifnot(is.matrix(C), nrow(C) == P, ncol(C) == P)
  
  # Beta handling
  if (is.null(beta_mat)) {
    beta_mat <- matrix(0.16, nrow = n_times, ncol = P)
  } else if (length(beta_mat) == 1) {
    beta_mat <- matrix(as.numeric(beta_mat), nrow = n_times, ncol = P)
  } else if (is.vector(beta_mat) && length(beta_mat) == n_times) {
    beta_mat <- matrix(beta_mat, nrow = n_times, ncol = P)
  } else {
    stopifnot(is.matrix(beta_mat), nrow(beta_mat) == n_times, ncol(beta_mat) == P)
  }
  if (any(!is.finite(beta_mat)) || any(beta_mat < 0))
    stop("beta_mat must be finite and non-negative.")
  
  # epsilon/alpha shape
  if (length(epsilon) == 1) epsilon <- rep(epsilon, P)
  if (!is.null(alpha) && length(alpha) == 1) alpha <- rep(alpha, P)
  stopifnot(length(epsilon) == P)
  if (!is.null(alpha)) stopifnot(length(alpha) == P, all(alpha >= 0 & alpha <= 1))
  
  pr <- function(x) pmax(pmin(x, 1), 0)
  
  # States (counts): [time, sims, pops]
  S <- array(0L, dim = c(n_times, n_sims, P))
  I <- array(0L, dim = c(n_times, n_sims, P))
  R <- array(0L, dim = c(n_times, n_sims, P))
  cases <- array(0L, dim = c(n_times, n_sims, P))
  
  # Init t=1
  for (p in seq_len(P)) {
    S[1, , p] <- as.integer(pop_vec[p] - I_init[p])
    I[1, , p] <- as.integer(I_init[p])
    R[1, , p] <- 0L
  }
  
  for (t in 2:n_times) {
    # Prevalence by source pop at t-1: [sims x P]
    I_slice <- I[t - 1, , ]                 # [n_sims x P]
    prev <- sweep(I_slice, 2, pop_vec, "/") # divide columns by N_p
    
    # Receiver-weighted prevalence per sim: [sims x P]
    # C rows=receivers, cols=sources → prev %*% t(C) gives receivers in columns
    M <- prev %*% t(C)
    
    beta_row <- beta_mat[t - 1, ]           # length P
    lambda_mat <- sweep(M, 2, beta_row, `*`)    # multiply each receiver col by its beta
    lambda_mat <- sweep(lambda_mat, 2, epsilon, `+`)  # add epsilon_p
    
    if (!stochastic) {
      new_inf <- matrix(0L, nrow = n_sims, ncol = P)
      new_rec <- matrix(0L, nrow = n_sims, ncol = P)
      loss_im <- matrix(0L, nrow = n_sims, ncol = P)
      for (p in seq_len(P)) {
        pinf <- pr(lambda_mat[, p])
        new_inf[, p] <- pmin(S[t-1, , p], round(pinf * S[t-1, , p]))
        new_rec[, p] <- round(gamma * I[t-1, , p])
        loss_im[, p] <- round(omega * R[t-1, , p])
      }
    } else {
      p_inf <- 1 - exp(-pr(lambda_mat))     # [sims x P]
      new_inf <- matrix(0L, nrow = n_sims, ncol = P)
      new_rec <- matrix(0L, nrow = n_sims, ncol = P)
      loss_im <- matrix(0L, nrow = n_sims, ncol = P)
      for (p in seq_len(P)) {
        new_inf[, p] <- stats::rbinom(n_sims, size = S[t-1, , p], prob = pr(p_inf[, p]))
        new_rec[, p] <- stats::rbinom(n_sims, size = I[t-1, , p], prob = pr(gamma))
        loss_im[, p] <- stats::rbinom(n_sims, size = R[t-1, , p], prob = pr(omega))
      }
    }
    
    # Update & record
    for (p in seq_len(P)) {
      S[t, , p] <- S[t-1, , p] - new_inf[, p] + loss_im[, p]
      I[t, , p] <- I[t-1, , p] + new_inf[, p] - new_rec[, p]
      R[t, , p] <- R[t-1, , p] - loss_im[, p] + new_rec[, p]
      
      cases[t, , p] <- if (is.null(alpha)) new_inf[, p]
      else stats::rbinom(n_sims, size = new_inf[, p], prob = pr(alpha[p]))
      
      S[t, , p] <- pmax(S[t, , p], 0L)
      I[t, , p] <- pmax(I[t, , p], 0L)
      R[t, , p] <- pmax(R[t, , p], 0L)
    }
  }
  
  # Proportions [time, sims, pops, state]
  props <- array(NA_real_, dim = c(n_times, n_sims, P, 3),
                 dimnames = list(NULL, NULL, NULL, c("S","I","R")))
  for (p in seq_len(P)) {
    props[, , p, "S"] <- S[, , p] / pop_vec[p]
    props[, , p, "I"] <- I[, , p] / pop_vec[p]
    props[, , p, "R"] <- R[, , p] / pop_vec[p]
  }
  
  list(
    time        = seq_len(n_times),
    proportions = props,
    cases       = cases,
    params      = list(
      n_times=n_times, pop_vec=pop_vec, I_init=I_init,
      beta=beta_mat, gamma=gamma, omega=omega,
      C=C, epsilon=epsilon, alpha=alpha,
      n_sims=n_sims, stochastic=stochastic, seed=seed
    )
  )
}
