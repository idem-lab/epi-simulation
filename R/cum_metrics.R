# Project: Kids Research Institute — SIRS modelling
# File: R/cum_metrics.R
# Purpose: Cumulative incidence & attack rate with robust shape handling

# safe-null coalescer
`%||%` <- function(x, y) if (is.null(x)) y else x

# Infer per-group populations (vector) or single scalar
.get_pop_vec <- function(sim) {
  # Common places we stored pop info
  pv  <- sim$params$pop_vec %||% sim$pop_vec
  p1  <- sim$params$pop     %||% sim$pop
  if (!is.null(pv)) return(as.numeric(pv))
  if (!is.null(p1)) return(as.numeric(p1))
  1
}

# ---- CUMULATIVE INCIDENCE ----------------------------------------------------

#' cumulative_incidence
#' @title Cumulative incidence over time
#' @description
#' Returns cumulative sums of daily infections (counts), handling shapes:
#' * deterministic single-pop: `incidence` vector
#' * stochastic single-pop: `cases` matrix [time × sims]
#' * multi-pop deterministic: `incidence` matrix [time × groups]
#' * multi-pop stochastic: `cases` **array** [time × sims × groups]
#' If a tidy helper `to_tidy()` exists and the structure is unusual, it will be
#' used as a fallback.
#' @param sim A simulation output.
#' @return `data.frame` with columns: `time`, `cum_incidence`, and `sim`/`group`
#'   when available.
#' @examples
#' # df <- cumulative_incidence(sim)
cumulative_incidence <- function(sim) {
  # ---- 1) Use CASES if present (preferred for reported flows) ----------------
  if (!is.null(sim$cases)) {
    cs <- sim$cases
    
    # [time × sims]
    if (is.matrix(cs)) {
      out <- lapply(seq_len(ncol(cs)), function(j) {
        data.frame(
          time = sim$time,
          sim  = j,
          cum_incidence = cumsum(cs[, j])
        )
      })
      return(do.call(rbind, out))
    }
    
    # [time × sims × groups]
    if (length(dim(cs)) == 3L) {
      TT <- dim(cs)[1]; M <- dim(cs)[2]; G <- dim(cs)[3]
      dn <- dimnames(cs)
      sims   <- dn[[2]] %||% as.character(seq_len(M))
      groups <- dn[[3]] %||% as.character(seq_len(G))
      
      rows <- vector("list", M * G); k <- 0L
      for (g in seq_len(G)) {
        for (m in seq_len(M)) {
          k <- k + 1L
          rows[[k]] <- data.frame(
            time = sim$time %||% seq_len(TT),
            sim  = sims[m],
            group = groups[g],
            cum_incidence = cumsum(cs[, m, g])
          )
        }
      }
      return(do.call(rbind, rows))
    }
  }
  
  # ---- 2) Otherwise use INCIDENCE if present ---------------------------------
  if (!is.null(sim$incidence)) {
    inc <- sim$incidence
    
    # vector [time]
    if (is.vector(inc)) {
      return(data.frame(
        time = sim$time,
        cum_incidence = cumsum(inc)
      ))
    }
    
    # [time × groups]
    if (is.matrix(inc)) {
      G <- ncol(inc)
      groups <- colnames(inc) %||% as.character(seq_len(G))
      out <- lapply(seq_len(G), function(g) {
        data.frame(
          time = sim$time,
          group = groups[g],
          cum_incidence = cumsum(inc[, g])
        )
      })
      return(do.call(rbind, out))
    }
    
    # [time × sims × groups]
    if (length(dim(inc)) == 3L) {
      TT <- dim(inc)[1]; M <- dim(inc)[2]; G <- dim(inc)[3]
      dn <- dimnames(inc)
      sims   <- dn[[2]] %||% as.character(seq_len(M))
      groups <- dn[[3]] %||% as.character(seq_len(G))
      
      rows <- vector("list", M * G); k <- 0L
      for (g in seq_len(G)) {
        for (m in seq_len(M)) {
          k <- k + 1L
          rows[[k]] <- data.frame(
            time = sim$time %||% seq_len(TT),
            sim  = sims[m],
            group = groups[g],
            cum_incidence = cumsum(inc[, m, g])
          )
        }
      }
      return(do.call(rbind, rows))
    }
  }
  
  # ---- 3) Fallback: try to_tidy() --------------------------------------------
  if (exists("to_tidy") && is.function(to_tidy)) {
    td <- try(to_tidy(sim), silent = TRUE)
    if (!inherits(td, "try-error")) {
      # Prefer explicit daily infections if available
      daily_col <- intersect(c("incidence", "cases"), names(td))
      if (length(daily_col)) {
        val <- daily_col[1]
        group_cols <- intersect(c("group", "sim"), names(td))
        td <- td[td[[val]] >= 0 & is.finite(td[[val]]), , drop = FALSE]
        td <- td[order(td$time), , drop = FALSE]
        agg <- aggregate(td[[val]],
                         by = c(list(td$time), lapply(group_cols, function(nm) td[[nm]])),
                         FUN = sum, na.rm = TRUE)
        names(agg) <- c("time", group_cols, "x")
        # cumulative by group_cols
        do_by <- if (length(group_cols)) group_cols else NULL
        if (length(do_by)) {
          res <- do.call(rbind, lapply(split(agg, agg[do_by]), function(df) {
            df$cum_incidence <- cumsum(df$x); df$x <- NULL; df
          }))
        } else {
          agg$cum_incidence <- cumsum(agg$x); agg$x <- NULL
          res <- agg
        }
        return(res)
      }
    }
  }
  
  stop("Unrecognized sim structure for cumulative_incidence().")
}

# ---- ATTACK RATE -------------------------------------------------------------

#' attack_rate
#' @title C