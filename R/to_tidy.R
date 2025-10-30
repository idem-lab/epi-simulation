#' Tidy a SIRS simulation into a long data frame
#'
#' @description
#' Convert heterogeneous SIRS simulation outputs (deterministic or stochastic,
#' single- or multi-population) into one **long** table with a consistent schema:
#' columns `time`, `group`, `sim`, `state`, `value`.
#'
#' @param sim A simulation result list (deterministic or stochastic; single- or multi-pop)
#'   produced by the package simulators. See **Supported input shapes**.
#'
#' @return A base `data.frame` in long (tidy) format with columns
#'   `time, group, sim, state, value`.
#'
#' @examples
#' # Deterministic single-pop
#' # df <- to_tidy(sim_det_single); head(df)
#'
#' # Deterministic multi-pop
#' # df <- to_tidy(sim_det_multi)
#'
#' # Stochastic single-pop
#' # df <- to_tidy(sim_stoch_single)
#'
#' # Stochastic multi-pop
#' # df <- to_tidy(sim_stoch_multi)
#'
#' @export


to_tidy <- function(sim) {
  pack <- function(t, g, s, S, I, R, inc) {
    rbind(
      data.frame(time = t, group = g, sim = s, state = "S",          value = S),
      data.frame(time = t, group = g, sim = s, state = "I",          value = I),
      data.frame(time = t, group = g, sim = s, state = "R",          value = R),
      data.frame(time = t, group = g, sim = s, state = "incidence",  value = inc)
    )
  }
  
  # ---------- Deterministic single-pop (vectors S/I/R/incidence) ----------
  if (!is.null(sim$S) && is.vector(sim$S) && is.null(sim$proportions)) {
    return(pack(sim$time, 1, 1, sim$S, sim$I, sim$R, sim$incidence))
  }
  
  # ---------- Deterministic multi-pop (matrices [T x P]) ----------
  if (!is.null(sim$S) && is.matrix(sim$S) && !is.null(sim$incidence) && is.matrix(sim$incidence) && is.null(sim$proportions)) {
    G <- ncol(sim$S)
    out <- lapply(seq_len(G), function(g) {
      pack(sim$time, g, 1, sim$S[, g], sim$I[, g], sim$R[, g], sim$incidence[, g])
    })
    return(do.call(rbind, out))
  }
  
  # ---------- Stochastic cases ----------
  if (!is.null(sim$proportions)) {
    prop <- sim$proportions
    d    <- dim(prop)
    dn   <- dimnames(prop)
    
    # helper: map state indices
    state_names <- if (length(dn) >= 3) dn[[length(dn)]] else NULL
    get_si <- function(nm) {
      if (is.null(state_names)) {
        # fallback assume order S,I,R
        c(S = 1L, I = 2L, R = 3L)[[nm]]
      } else {
        idx <- match(nm, state_names)
        if (is.na(idx)) stop("State '", nm, "' not found in proportions dimnames.")
        idx
      }
    }
    iS <- get_si("S"); iI <- get_si("I"); iR <- get_si("R")
    
    # Cases shape can be [T x sims] (single-pop) or [T x sims x P] (multi-pop)
    cases <- sim$cases
    
    if (length(d) == 3L) {
      # [T x sims x state]  -> single-pop
      Tn <- d[1]; M <- d[2]
      if (!is.null(dim(cases)) && !all(dim(cases)[1:2] == c(Tn, M))) {
        stop("cases shape doesn't match proportions for single-pop.")
      }
      out <- lapply(seq_len(M), function(m) {
        pack(sim$time, 1, m,
             S   = prop[, m, iS],
             I   = prop[, m, iI],
             R   = prop[, m, iR],
             inc = if (is.null(dim(cases))) cases else cases[, m])
      })
      return(do.call(rbind, out))
    }
    
    if (length(d) == 4L) {
      # [T x sims x P x state]  -> multi-pop
      Tn <- d[1]; M <- d[2]; P <- d[3]
      if (is.null(dim(cases)) || !all(dim(cases)[1:3] == c(Tn, M, P))) {
        stop("cases shape must be [T x sims x P] for multi-pop stochastic.")
      }
      out <- lapply(seq_len(P), function(g) {
        lapply(seq_len(M), function(m) {
          pack(sim$time, g, m,
               S   = prop[, m, g, iS],
               I   = prop[, m, g, iI],
               R   = prop[, m, g, iR],
               inc = cases[, m, g])
        })
      })
      return(do.call(rbind, unlist(out, recursive = FALSE)))
    }
    
    stop("Unsupported proportions shape: expected 3 or 4 dims.")
  }
  
  stop("Unrecognized sim structure for to_tidy().")
}









