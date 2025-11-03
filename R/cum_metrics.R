# =============================================================================
# Project: Kids Research Institute — SIRS modelling
# File   : R/cum_metrics.R
# Purpose: Cumulative incidence and attack-rate helpers
# Inputs : sim (deterministic | stochastic | multi-pop; single or multi-run)
# Outputs:
#   - cumulative_incidence(sim) -> data.frame with time and cum_* plus
#     optional columns: sim (simulation index) and/or group (group index)
#   - attack_rate(sim)          -> scalar / vector / matrix (see below)
# Notes : Pure base R (no dplyr). ASCII only.
# =============================================================================

# Small helper: x %||% y  -> return x unless it's NULL
`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal: extract population vector per group (or scalar for single-pop)
.get_pop_vec <- function(sim) {
  # Prefer explicit vector if present
  pv <- sim$params$pop_vec %||% sim$pop_vec
  if (!is.null(pv)) return(as.numeric(pv))

  # Fall back to scalar pop in params or top-level
  p  <- sim$params$pop %||% sim$pop
  if (!is.null(p)) return(as.numeric(p))

  # Last resort: 1 (prevents division by zero)
  1
}

# -----------------------------------------------------------------------------
# cumulative_incidence()
# -----------------------------------------------------------------------------
# Returns cumulative sums of daily infections (preferred: 'incidence').
# If 'incidence' is absent, falls back to 'cases' (reported).
# Shapes handled:
# - deterministic single-pop: vector [time]
# - deterministic multi-pop : matrix [time x groups]
# - stochastic single-pop   : matrix [time x sims]
# - stochastic multi-pop    : array  [time x sims x groups]
#
# Output is always a data.frame with columns:
# - time
# - cum_incidence  (or cum_cases if incidence absent)
# - optional: sim (for stochastic), group (for multi-pop)
# -----------------------------------------------------------------------------
cumulative_incidence <- function(sim) {
  use_field <- if (!is.null(sim$incidence)) "incidence" else if (!is.null(sim$cases)) "cases" else NULL
  if (is.null(use_field)) stop("No 'incidence' or 'cases' found in sim for cumulative_incidence().")

  x <- sim[[use_field]]
  nm <- if (use_field == "incidence") "cum_incidence" else "cum_cases"

  # 1) Vector [time]
  if (is.vector(x) && is.null(dim(x))) {
    return(data.frame(
      time = sim$time,
      ..tmp.. = cumsum(as.numeric(x))
    , check.names = FALSE, fix.empty.names = FALSE,
    row.names = NULL, stringsAsFactors = FALSE)[, c("time", "..tmp.."),
      drop = FALSE
    ] -> out) # little trick to avoid partial matching
    names(out) <- c("time", nm)
    return(out)
  }

  # 2) Matrix [time x K]
  if (is.matrix(x)) {
    # Determine whether K is sims or groups by checking other fields
    K <- ncol(x)
    # Heuristic: if sim$sims exists or sim$params$n_sims matches K -> sims
    n_sims <- sim$params$n_sims %||% sim$n_sims
    is_sims <- !is.null(n_sims) && as.numeric(n_sims) == K

    # Build DF without mixing columns
    if (is_sims) {
      # single-pop stochastic: sims are columns
      out_list <- vector("list", K)
      for (j in seq_len(K)) {
        out_list[[j]] <- data.frame(
          time = sim$time,
          sim  = j,
          ..tmp.. = cumsum(as.numeric(x[, j]))
        , check.names = FALSE, fix.empty.names = FALSE,
        row.names = NULL, stringsAsFactors = FALSE)
      }
      out <- do.call(rbind, out_list)
      names(out)[names(out) == "..tmp.."] <- nm
      rownames(out) <- NULL
      return(out)
    } else {
      # deterministic multi-pop: groups are columns
      out_list <- vector("list", K)
      for (g in seq_len(K)) {
        out_list[[g]] <- data.frame(
          time  = sim$time,
          group = g,
          ..tmp.. = cumsum(as.numeric(x[, g]))
        , check.names = FALSE, fix.empty.names = FALSE,
        row.names = NULL, stringsAsFactors = FALSE)
      }
      out <- do.call(rbind, out_list)
      names(out)[names(out) == "..tmp.."] <- nm
      rownames(out) <- NULL
      return(out)
    }
  }

  # 3) Array [time x sims x groups]
  if (length(dim(x)) == 3L) {
    TT <- dim(x)[1]; M <- dim(x)[2]; G <- dim(x)[3]
    out_list <- vector("list", M * G)
    idx <- 0L
    for (j in seq_len(M)) {
      for (g in seq_len(G)) {
        idx <- idx + 1L
        out_list[[idx]] <- data.frame(
          time  = sim$time,
          sim   = j,
          group = g,
          ..tmp.. = cumsum(as.numeric(x[, j, g]))
        , check.names = FALSE, fix.empty.names = FALSE,
        row.names = NULL, stringsAsFactors = FALSE)
      }
    }
    out <- do.call(rbind, out_list)
    names(out)[names(out) == "..tmp.."] <- nm
    rownames(out) <- NULL
    return(out)
  }

  stop("Unrecognized data shape in cumulative_incidence().")
}

# -----------------------------------------------------------------------------
# attack_rate()
# -----------------------------------------------------------------------------
# Cumulative attack rate over the simulation horizon:
# - Prefer 'cases' (reported) if present; else 'incidence' (true infections).
# - Divides by population (scalar for single-pop or vector per group).
# - Shapes:
#   * single-pop deterministic (vector) -> scalar
#   * single-pop stochastic [time x sims] -> length M vector
#   * multi-pop deterministic [time x groups] -> length G vector
#   * multi-pop stochastic [time x sims x groups] -> matrix [M x G]
# - Fallback (warning): final R proportion(s) if no flows available.
# -----------------------------------------------------------------------------
attack_rate <- function(sim) {
  pop_vec <- .get_pop_vec(sim)

  # ---------- Prefer CASES ----------
  if (!is.null(sim$cases)) {
    cs <- sim$cases

    # [time] -> scalar
    if (is.vector(cs) && is.null(dim(cs))) {
      pop <- if (length(pop_vec) == 1L) pop_vec else pop_vec[1L]
      return(sum(as.numeric(cs), na.rm = TRUE) / pop)
    }

    # [time x sims] -> vector length sims
    if (is.matrix(cs)) {
      pop <- if (length(pop_vec) == 1L) pop_vec else pop_vec[1L]
      return(colSums(cs, na.rm = TRUE) / pop)
    }

    # [time x sims x groups] -> matrix [sims x groups]
    if (length(dim(cs)) == 3L) {
      sums <- apply(cs, c(2, 3), sum, na.rm = TRUE)  # [M x G]
      G <- dim(cs)[3]
      pop_g <- if (length(pop_vec) == 1L) rep(pop_vec, G) else as.numeric(pop_vec)
      return(sweep(sums, 2, pop_g, "/"))
    }
  }

  # ---------- Then INCIDENCE ----------
  if (!is.null(sim$incidence)) {
    inc <- sim$incidence

    # [time] -> scalar
    if (is.vector(inc) && is.null(dim(inc))) {
      pop <- if (length(pop_vec) == 1L) pop_vec else pop_vec[1L]
      return(sum(as.numeric(inc), na.rm = TRUE) / pop)
    }

    # [time x groups] -> vector length groups
    if (is.matrix(inc)) {
      G <- ncol(inc)
      pop_g <- if (length(pop_vec) == 1L) rep(pop_vec, G) else as.numeric(pop_vec)
      return(colSums(inc, na.rm = TRUE) / pop_g)
    }

    # [time x sims x groups] -> matrix [sims x groups]
    if (length(dim(inc)) == 3L) {
      sums <- apply(inc, c(2, 3), sum, na.rm = TRUE)  # [M x G]
      G <- dim(inc)[3]
      pop_g <- if (length(pop_vec) == 1L) rep(pop_vec, G) else as.numeric(pop_vec)
      return(sweep(sums, 2, pop_g, "/"))
    }
  }

  # ---------- Fallback: final R proportion(s) ----------
  warning("No cases/incidence available; returning final R proportion as proxy (may understate true AR in SIRS).")

  # single-pop deterministic: R is vector
  if (!is.null(sim$R) && is.vector(sim$R) && is.null(sim$proportions)) {
    return(tail(as.numeric(sim$R), 1))
  }

  # single-pop stochastic: proportions [time x sims x state]
  if (!is.null(sim$proportions)) {
    M <- dim(sim$proportions)[2]
    out <- numeric(M)
    for (j in seq_len(M)) {
      out[j] <- tail(as.numeric(sim$proportions[, j, "R"]), 1)
    }
    return(out)
  }

  # multi-pop deterministic: R [time x groups]
  if (!is.null(sim$R) && is.matrix(sim$R)) {
    G <- ncol(sim$R)
    out <- numeric(G)
    for (g in seq_len(G)) {
      out[g] <- tail(as.numeric(sim$R[, g]), 1)
    }
    return(out)
  }

  stop("Unrecognized sim structure for attack_rate().")
}
