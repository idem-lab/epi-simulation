#' @title Summarize key outbreak metrics (single- or multi-pop)
#' @description
#' Extract headline diagnostics from SIRS simulations:
#' - peak infected **proportion** and its day,
#' - peak daily **incidence (count)** and its day,
#' - final recovered **proportion**.
#'
#' Supports:
#'  (A) Deterministic single-pop (vectors S/I/R, incidence)
#'  (A') Deterministic multi-pop (matrices [time x P] S/I/R, incidence)
#'  (B) Stochastic single-pop (proportions [time, sim, state], cases [time, sim])
#'  (B') Stochastic multi-pop (proportions [time, sim, P, state], cases [time, sim, P])
#'
#' @param sim   Simulation result list.
#' @param level One of "total" (aggregate across pops) or "per_pop" (one row per pop).
#'              Only relevant for multi-pop inputs. Default "total".
#'
#' @return A data.frame with columns:
#'   sims, (pop if per_pop), peak_I, peak_I_day, peak_incidence,
#'   peak_incidence_day, final_R.
#' @export
summarize_sim <- function(sim, level = c("total","per_pop")) {
  level <- match.arg(level)
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # helper to compute metrics for a single I-series (prop) and incidence-series (count)
  .summ <- function(I_vec, inc_vec, R_last, sim_id = 1L, pop_id = NA_integer_) {
    I_vec[!is.finite(I_vec)]   <- NA_real_
    inc_vec[!is.finite(inc_vec)] <- NA_real_
    data.frame(
      sims = sim_id,
      pop  = pop_id,
      peak_I = max(I_vec, na.rm = TRUE),
      peak_I_day = which.max(I_vec),
      peak_incidence = max(inc_vec, na.rm = TRUE),
      peak_incidence_day = which.max(inc_vec),
      final_R = R_last,
      row.names = NULL
    )
  }
  
  # ---------------- Case A: deterministic single-pop (vectors) ----------------
  if (!is.null(sim$S) && is.vector(sim$I) && is.null(sim$proportions)) {
    return(.summ(sim$I, sim$incidence, tail(sim$R, 1)))
  }
  
  # -------------- Case A': deterministic multi-pop (matrices [T x P]) --------
  if (!is.null(sim$S) && is.matrix(sim$I) && is.null(sim$proportions)) {
    I_mat   <- sim$I            # proportions per pop
    R_mat   <- sim$R            # proportions per pop
    inc_mat <- sim$incidence    # counts per pop (assumed)
    T <- nrow(I_mat); P <- ncol(I_mat)
    
    # population sizes (if present) to do weighted totals; else equal weights
    pop_vec <- sim$params$pop_vec %||% rep(1, P)
    pop_tot <- sum(pop_vec)
    
    if (level == "per_pop") {
      out <- lapply(seq_len(P), function(p)
        .summ(I_mat[, p], inc_mat[, p], R_mat[T, p], sim_id = 1L, pop_id = p))
      return(do.call(rbind, out))
    } else {
      # weighted proportion totals, summed counts
      I_tot   <- as.vector(I_mat %*% (pop_vec / pop_tot))
      R_tot   <- as.vector(R_mat %*% (pop_vec / pop_tot))
      inc_tot <- rowSums(inc_mat)
      return(.summ(I_tot, inc_tot, R_tot[T], sim_id = 1L, pop_id = NA_integer_))
    }
  }
  
  # --------- Case B: stochastic single-pop (proportions [T, S, state]) -------
  if (!is.null(sim$proportions) && length(dim(sim$proportions)) == 3) {
    I_arr  <- sim$proportions[, , "I", drop = FALSE]  # [T, S, 1]
    R_arr  <- sim$proportions[, , "R", drop = FALSE]  # [T, S, 1]
    inc_m  <- sim$cases                                # [T, S]
    T <- dim(I_arr)[1]; S <- dim(I_arr)[2]
    out <- lapply(seq_len(S), function(j)
      .summ(I_arr[, j, 1], inc_m[, j], R_arr[T, j, 1], sim_id = j, pop_id = NA_integer_))
    return(do.call(rbind, out))
  }
  
  # ---- Case B': stochastic multi-pop (proportions [T, S, P, state]) ---------
  if (!is.null(sim$proportions) && length(dim(sim$proportions)) == 4) {
    # dims: [time, sim, pop, state]
    dims <- dim(sim$proportions)
    T <- dims[1]; S <- dims[2]; P <- dims[3]
    
    I_arr <- sim$proportions[, , , "I", drop = FALSE]  # [T, S, P, 1]
    R_arr <- sim$proportions[, , , "R", drop = FALSE]  # [T, S, P, 1]
    inc_a <- sim$cases                                 # [T, S, P] counts
    
    pop_vec <- sim$params$pop_vec %||% rep(1, P)
    pop_tot <- sum(pop_vec)
    w <- pop_vec / pop_tot
    
    if (level == "per_pop") {
      out <- lapply(seq_len(S), function(j) {
        lapply(seq_len(P), function(p) {
          Ijp   <- I_arr[, j, p, 1]
          Rjp   <- R_arr[, j, p, 1]
          incjp <- inc_a[, j, p]
          .summ(Ijp, incjp, Rjp[T], sim_id = j, pop_id = p)
        })
      })
      return(do.call(rbind, unlist(out, recursive = FALSE)))
    } else {
      # totals across pops: weighted proportions, summed counts
      out <- lapply(seq_len(S), function(j) {
        Ij   <- as.vector(I_arr[, j, , 1, drop = TRUE] %*% w)  # [T]
        Rj   <- as.vector(R_arr[, j, , 1, drop = TRUE] %*% w)  # [T]
        incj <- rowSums(inc_a[, j, , drop = TRUE])             # [T]
        .summ(Ij, incj, Rj[T], sim_id = j, pop_id = NA_integer_)
      })
      return(do.call(rbind, out))
    }
  }
  
  stop("summarize_sim: Input shape not recognized. Please check your sim object.")
}
