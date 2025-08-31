#' cumulative_incidence
#' @title Cumulative incidence over time
#' @description
#' Returns cumulative sums of daily infections (counts), handling each sim shape:
#' - deterministic single-pop (`incidence` vector)
#' - stochastic (`cases` matrix [time x sims])
#' - multi-pop (`incidence` matrix [time x groups])
#' @param sim A simulation output.
#' @return data.frame with `time` and `cum_incidence` (+ `sim`/`group` if present).
#' @examples
#' # ci <- cumulative_incidence(sim); head(ci)
cumulative_incidence <- function(sim) {
  # Deterministic single-pop
  if (!is.null(sim$incidence) && is.vector(sim$incidence) && is.null(sim$cases)) {
    return(data.frame(time = sim$time, cum_incidence = cumsum(sim$incidence)))
  }
  # Stochastic multi-run
  if (!is.null(sim$cases)) {
    ci <- apply(sim$cases, 2, cumsum)  # [time x sims]
    out <- lapply(seq_len(ncol(ci)), function(j) {
      data.frame(time = sim$time, sim = j, cum_incidence = ci[, j])
    })
    return(do.call(rbind, out))
  }
  # Multi-pop deterministic
  if (!is.null(sim$incidence) && is.matrix(sim$incidence)) {
    ci <- apply(sim$incidence, 2, cumsum)  # [time x groups]
    out <- lapply(seq_len(ncol(ci)), function(g) {
      data.frame(time = sim$time, group = g, cum_incidence = ci[, g])
    })
    return(do.call(rbind, out))
  }
  stop("Unrecognized sim structure for cumulative_incidence().")
}

#' attack_rate
#' @title Final attack rate (approx. proportion infected at the end)
#' @description
#' Uses the final recovered proportion R(T) as an attack-rate proxy.
#' Returns a scalar (single-pop), a vector over sims (stochastic), or a vector
#' over groups (multi-pop).
#' @param sim A simulation output.
#' @return Numeric scalar or vector of final R proportions.
#' @examples
#' # ar <- attack_rate(sim)
attack_rate <- function(sim) {
  # Deterministic single-pop
  if (!is.null(sim$R) && is.vector(sim$R) && is.null(sim$proportions)) {
    return(tail(sim$R, 1))
  }
  # Stochastic multi-run
  if (!is.null(sim$proportions)) {
    M <- dim(sim$proportions)[2]
    return(sapply(seq_len(M), function(j) tail(sim$proportions[, j, "R"], 1)))
  }
  # Multi-pop deterministic
  if (!is.null(sim$R) && is.matrix(sim$R)) {
    G <- ncol(sim$R)
    return(sapply(seq_len(G), function(g) tail(sim$R[, g], 1)))
  }
  stop("Unrecognized sim structure for attack_rate().")
}
