#' @title Plot multi-population SIRS outputs (deterministic + stochastic)
#'
#' @description
#' Visualise SIRS simulation results for **multiple populations** from either
#' deterministic runs (`[T × P]` matrices) or stochastic runs
#' (`proportions [T × sims × P × {S,I,R}]`, `cases [T × sims × P]`).
#'
#' Panels supported:
#' - **S / I / R panels**
#'   - `group_style = "facet"` → facet by population; fixed S/I/R colours.
#'   - `group_style = "combined"` → one panel; colours map to populations.
#' - **Incidence panel** (always combined): mean line per population with optional
#'   quantile ribbon across simulations.
#'
#' @param sim See description (det matrices or stoch arrays) with `time`.
#' @param probs Length-2 numeric for ribbon quantiles, or single `p` in (0,1).
#' @param per_million Logical; if `TRUE`, scale incidence by population (per million).
#' @param which One of `c("S","I","R","incidence")`.
#' @param group_style `"facet"` or `"combined"` for S/I/R panels.
#' @param pop_cols Optional vector of population colours.
#' @param show_bands Logical; draw quantile ribbons (stochastic only).
#'
#' @export
plot_multi <- function(
    sim,
    probs       = c(0.025, 0.975),
    per_million = FALSE,
    which       = c("S", "I", "R", "incidence"),
    group_style = c("facet", "combined"),
    pop_cols    = NULL,
    show_bands  = TRUE
){
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_multi() needs {ggplot2}.")
  which       <- match.arg(which)
  group_style <- match.arg(group_style)
  
  # ------- helpers-------
  `%||%` <- function(x, y) if (is.null(x)) y else x
  .not_blank_df <- function(df) !is.null(df) && is.data.frame(df) && nrow(df) > 0
  
  .pop_cols <- function(P) {
    c("#3B528B", "#5DC863", "#FDE725", "#21918C", "#440154", "#F8961E", "#58A4B0")[seq_len(min(P,7))]
  }
  STATE_COLS <- c(S = "#0072B2", I = "#CC79A7", R = "#009E73")
  .alpha <- function(cols, a = 0.25) {
    if (length(cols) == 1L) grDevices::adjustcolor(cols, a) else vapply(cols, grDevices::adjustcolor, "", alpha.f = a)
  }
  .normalize_pop_names <- function(P, pop_names = NULL){
    if (is.null(pop_names)) return(paste0("Pop ", seq_len(P)))
    if (length(pop_names) == P) return(pop_names)
    warning(sprintf("Length of pop_names (%d) != P (%d); using defaults.", length(pop_names), P))
    paste0("Pop ", seq_len(P))
  }
  .norm_probs <- function(probs){
    if (missing(probs) || is.null(probs)) return(c(0.025, 0.975))
    stopifnot(is.numeric(probs), all(is.finite(probs)))
    if (length(probs) == 1L) {
      p <- as.numeric(probs); stopifnot(p > 0, p < 1)
      return(c(1 - p, p))
    }
    stopifnot(length(probs) >= 2)
    sort(probs[1:2])
  }
  .to_mat_lenT <- function(x, T){
    if (is.null(x)) return(NULL)
    if (is.matrix(x)){ stopifnot(nrow(x) == T); return(x) }
    x <- as.numeric(x); stopifnot(length(x) == T); matrix(x, ncol = 1)
  }
  .ensure_4d_proportions <- function(prop, T){
    if (is.null(prop)) return(NULL)
    d <- dim(prop)
    if (length(d) == 4L){
      stopifnot(d[1] == T)
      dn <- dimnames(prop)
      if (is.null(dn) || is.null(dn[[4]])){
        k <- d[4]; labs <- c("S","I","R","X","Y","Z")[seq_len(k)]
        dimnames(prop) <- if (is.null(dn)) list(NULL,NULL,NULL,labs) else { dn[[4]] <- labs; dn }
      }
      return(prop)
    }
    if (length(d) == 3L){
      stopifnot(d[1] == T)
      k <- d[3]; labs <- dimnames(prop); st <- if (!is.null(labs)) labs[[3]] else NULL
      if (is.null(st)) st <- c("S","I","R")[seq_len(k)]
      array(prop, dim = c(d[1], d[2], 1L, k), dimnames = list(NULL,NULL,NULL,st))
    } else stop("`proportions` must be 3D [T x sims x state] or 4D [T x sims x pop x state].")
  }
  .ensure_3d_cases <- function(cases, T){
    if (is.null(cases)) return(NULL)
    d <- dim(cases)
    if (length(d) == 3L){ stopifnot(d[1] == T); return(cases) }
    if (length(d) == 2L){ stopifnot(d[1] == T); return(array(cases, dim = c(d[1], d[2], 1L))) }
    stop("`cases` must be 2D [T x sims] or 3D [T x sims x pop].")
  }
  .coerce_sim_for_dashboard <- function(sim){
    stopifnot(!is.null(sim$time))
    T <- length(sim$time)
    S <- .to_mat_lenT(sim$S,T); I <- .to_mat_lenT(sim$I,T); R <- .to_mat_lenT(sim$R,T); Inc <- .to_mat_lenT(sim$incidence,T)
    prop <- .ensure_4d_proportions(sim$proportions,T)
    cases<- .ensure_3d_cases(sim$cases,T)
    P_det <- suppressWarnings(max(c(ncol(S), ncol(I), ncol(R), ncol(Inc)), na.rm = TRUE)); if (!is.finite(P_det)) P_det <- NA_integer_
    P_sto <- if (!is.null(prop)) dim(prop)[3] else NA_integer_
    if (!is.na(P_sto) && !is.na(P_det) && P_det != P_sto){
      if (P_det == 1 && P_sto > 1){
        if (!is.null(S))   S   <- matrix(S[,1],   T, P_sto)
        if (!is.null(I))   I   <- matrix(I[,1],   T, P_sto)
        if (!is.null(R))   R   <- matrix(R[,1],   T, P_sto)
        if (!is.null(Inc)) Inc <- matrix(Inc[,1], T, P_sto)
        P_det <- P_sto
      } else if (P_sto == 1 && P_det > 1){
        prop  <- array(prop,  dim = c(T, dim(prop)[2],  P_det, dim(prop)[4])); for (p in seq_len(P_det)) prop[,,p,] <- prop[,,1,]
        cases <- array(cases, dim = c(T, dim(cases)[2], P_det));              for (p in seq_len(P_det)) cases[,,p]   <- cases[,,1]
        P_sto <- P_det
      } else stop(sprintf("Population mismatch (det P=%s vs stoch P=%s).", P_det, P_sto))
    }
    P <- if (!is.na(P_det)) P_det else if (!is.na(P_sto)) P_sto else 1L
    params <- sim$params %||% list()
    params$pop_vec   <- if (is.null(params$pop_vec) || length(params$pop_vec) != P) rep(1, P) else params$pop_vec
    params$pop_names <- .normalize_pop_names(P, params$pop_names)
    params$beta      <- params$beta %||% NULL
    list(time=sim$time, S=S, I=I, R=R, incidence=Inc, proportions=prop, cases=cases, params=params)
  }
  .melt_mat <- function(mat, time, colname, pop_names=NULL){
    if (is.null(mat)) return(NULL)
    stopifnot(is.matrix(mat), nrow(mat) == length(time))
    T <- nrow(mat); P <- ncol(mat)
    df <- data.frame(
      time = rep(time, times = P),
      pop  = factor(rep(seq_len(P), each = T), labels = pop_names %||% paste0("Pop ", seq_len(P))),
      tmp  = as.vector(mat),
      check.names = FALSE
    )
    names(df)[3] <- colname   # <- create a column literally named "S", "I", "R", or "cases"
    df
  }
  .melt_band <- function(arr, time, probs = c(0.025, 0.975), per_million = FALSE, pop_size = NULL, pop_names = NULL){
    probs <- .norm_probs(probs)
    stopifnot(length(dim(arr)) == 3, dim(arr)[1] == length(time))
    T <- dim(arr)[1]; P <- dim(arr)[3]
    qL <- apply(arr, c(1,3), stats::quantile, probs[1], na.rm = TRUE)
    qU <- apply(arr, c(1,3), stats::quantile, probs[2], na.rm = TRUE)
    mu <- apply(arr, c(1,3), mean, na.rm = TRUE)
    if (per_million && !is.null(pop_size)) {
      scale <- matrix(pop_size/1e6, nrow = T, ncol = P, byrow = TRUE)
      mu <- mu/scale; qL <- qL/scale; qU <- qU/scale
    }
    data.frame(
      time = rep(time, times = P),
      pop  = factor(rep(seq_len(P), each = T), labels = pop_names %||% paste0("Pop ", seq_len(P))),
      mean = as.vector(mu),
      qL   = as.vector(qL),
      qU   = as.vector(qU),
      check.names = FALSE
    )
  }
  
  .build_SIR_panel <- function(state_name, det_df, band_df, ylab, P,
                               group_style = c("facet","combined"), cols = NULL, show_bands = TRUE){
    stopifnot(state_name %in% c("S","I","R"))
    group_style <- match.arg(group_style)
    state_col <- unname(STATE_COLS[[state_name]])
    single_pop <- (P == 1)
    
    # ---- combined branch ----
    if (group_style == "combined" && !single_pop){
      p <- ggplot2::ggplot()
      
      if (.not_blank_df(band_df)) {
        if (show_bands) {
          p <- p + ggplot2::geom_ribbon(
            data = band_df,
            ggplot2::aes(time, ymin = qL, ymax = qU, fill = pop),
            alpha = 0.25, inherit.aes = FALSE
          )
        }
        p <- p + ggplot2::geom_line(
          data = band_df, ggplot2::aes(time, mean, color = pop),
          linewidth = 0.9, inherit.aes = FALSE
        )
      }
      
      if (.not_blank_df(det_df)) {
        p <- p + ggplot2::geom_line(
          data = det_df, ggplot2::aes(time, .data[[state_name]], color = pop),
          linewidth = 0.9, inherit.aes = FALSE
        )
      }
      
      return(
        p +
          ggplot2::scale_color_manual(values = cols, drop = FALSE, name = "Population") +
          ggplot2::scale_fill_manual(values = stats::setNames(.alpha(cols, 0.25), names(cols)),
                                     drop = FALSE, name = "Population") +
          ggplot2::labs(title = paste0(state_name, "(t)"), x = "Day", y = ylab) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      )
    }
    
    # facet or single-pop branch
    p <- ggplot2::ggplot()
    if (.not_blank_df(band_df)) {
      if (show_bands) {
        p <- p + ggplot2::geom_ribbon(data = band_df, ggplot2::aes(time, ymin = qL, ymax = qU),
                                      inherit.aes = FALSE, fill = grDevices::adjustcolor(state_col, 0.25))
      }
      p <- p + ggplot2::geom_line(data = band_df, ggplot2::aes(time, mean, group = pop),
                                  inherit.aes = FALSE, linewidth = 0.9, color = state_col)
    }
    if (.not_blank_df(det_df)) {
      p <- p + ggplot2::geom_line(data = det_df, ggplot2::aes(time, .data[[state_name]]),
                                  inherit.aes = FALSE, linewidth = 0.9, color = state_col)
    }
    if (!single_pop && (.not_blank_df(band_df) || .not_blank_df(det_df))) {
      p <- p + ggplot2::facet_wrap(~pop, ncol = 2)
    }
    p + ggplot2::labs(title = paste0(state_name, "(t)"), x = "Day", y = ylab) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  # ------- main -------
  probs <- .norm_probs(probs)
  sim   <- .coerce_sim_for_dashboard(sim)
  time  <- sim$time
  
  has_det <- !is.null(sim$I) && is.matrix(sim$I)
  has_sto <- !is.null(sim$proportions) && length(dim(sim$proportions)) == 4L
  if (!has_det && !has_sto) stop("Unsupported `sim` shape: need deterministic matrices or stochastic arrays.")
  
  Tlen <- length(time)
  if (has_det) {
    for (nm in c("S","I","R","incidence")) {
      if (!is.null(sim[[nm]])) {
        if (!is.matrix(sim[[nm]])) stop(nm, " must be a [T × P] matrix.")
        if (nrow(sim[[nm]]) != Tlen) stop("Time length and ", nm, " rows differ (", Tlen, " vs ", nrow(sim[[nm]]), ").")
      }
    }
  } else {
    if (dim(sim$proportions)[1] != Tlen) stop("Time length and proportions first dim differ.")
    if (!is.null(sim$cases) && dim(sim$cases)[1] != Tlen) stop("Time length and cases first dim differ.")
  }
  
  P         <- if (has_det) ncol(sim$I) else dim(sim$proportions)[3]
  pop_names <- .normalize_pop_names(P, sim$params$pop_names)
  pops      <- sim$params$pop_vec %||% rep(1, P)
  
  # ----- palette -----
  .resolve_pop_cols <- function(P, pop_names, pop_cols_arg = NULL){
    # 1) explicit argument wins
    cols <- pop_cols_arg
    # 2) else global/project-level my_pop_cols if present
    if (is.null(cols)) cols <- get0("my_pop_cols", envir = parent.frame(), inherits = TRUE, ifnotfound = NULL)
    # 3) else fallback palette
    if (is.null(cols)) cols <- .pop_cols(P)
    cols <- unname(as.character(cols))
    if (length(cols) < P) cols <- rep(cols, length.out = P)
    if (length(cols) > P) cols <- cols[seq_len(P)]
    stats::setNames(cols, pop_names)
  }
  cols <- .resolve_pop_cols(P, pop_names, pop_cols_arg = pop_cols)
  
  inc_lab <- if (per_million) "Daily cases (per million)" else "Daily cases"
  s_lab <- "S proportion"; i_lab <- "I proportion"; r_lab <- "R proportion"
  
  # deterministic melts
  S_det <- I_det <- R_det <- Inc_det <- NULL
  if (has_det) {
    S_det <- .melt_mat(sim$S, time, "S", pop_names)
    I_det <- .melt_mat(sim$I, time, "I", pop_names)
    R_det <- .melt_mat(sim$R, time, "R", pop_names)
    inc_mat <- sim$incidence
    if (!is.null(inc_mat)) {
      if (per_million) inc_mat <- sweep(inc_mat, 2, pops/1e6, "/")
      Inc_det <- .melt_mat(inc_mat, time, "cases", pop_names)
    }
  }
  
  # stochastic bands
  I_band <- Inc_band <- S_band <- R_band <- NULL
  if (has_sto) {
    .get_state_3d <- function(prop4, state) {
      d <- dim(prop4); dn4 <- dimnames(prop4); stn <- if (is.null(dn4)) NULL else dn4[[4]]
      idx <- if (!is.null(stn) && state %in% stn) which(stn == state)[1] else {
        i <- switch(state, S = 1L, I = 2L, R = 3L, 1L); if (i > d[4]) 1L else i
      }
      array(prop4[,,,idx, drop = FALSE], dim = d[1:3])
    }
    S_arr  <- .get_state_3d(sim$proportions, "S")
    I_arr  <- .get_state_3d(sim$proportions, "I")
    R_arr  <- .get_state_3d(sim$proportions, "R")
    cases_a <- sim$cases
    if (length(dim(cases_a)) == 2L) cases_a <- array(cases_a, dim = c(dim(cases_a), 1L))
    
    S_band   <- .melt_band(S_arr,   time, probs, FALSE,   NULL, pop_names)
    I_band   <- .melt_band(I_arr,   time, probs, FALSE,   NULL, pop_names)
    R_band   <- .melt_band(R_arr,   time, probs, FALSE,   NULL, pop_names)
    Inc_band <- .melt_band(cases_a, time, probs, per_million, pops, pop_names)
  }
  
  # build panels
  plt_S <- .build_SIR_panel("S", S_det, S_band, s_lab, P,
                            group_style = group_style, cols = cols, show_bands = show_bands)
  plt_I <- .build_SIR_panel("I", I_det, I_band, i_lab, P,
                            group_style = group_style, cols = cols, show_bands = show_bands)
  plt_R <- .build_SIR_panel("R", R_det, R_band, r_lab, P,
                            group_style = group_style, cols = cols, show_bands = show_bands)
  
  # incidence
  single_pop <- (P == 1L)
  pplt <- ggplot2::ggplot()
  if (.not_blank_df(Inc_band)) {
    if (show_bands) {
      pplt <- pplt + ggplot2::geom_ribbon(
        ggplot2::aes(time, ymin = qL, ymax = qU, fill = pop),
        data = Inc_band, alpha = 0.25, show.legend = !single_pop
      )
    }
    pplt <- pplt + ggplot2::geom_line(
      ggplot2::aes(time, mean, color = pop),
      data = Inc_band, linewidth = 0.9, show.legend = !single_pop
    )
  }
  if (.not_blank_df(Inc_det)) {
    pplt <- pplt + ggplot2::geom_line(
      ggplot2::aes(time, cases, color = pop),
      data = Inc_det, linewidth = 0.9, show.legend = !single_pop
    )
  }
  plt_inc <- pplt +
    ggplot2::scale_color_manual(values = cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = stats::setNames(.alpha(cols, 0.25), names(cols)), drop = FALSE) +
    ggplot2::labs(title = "Daily Incidence", x = "Day", y = inc_lab, color = "Population", fill = "Population") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = if (single_pop) "none" else "right"
    )
  
  # return
  switch(which,
         S         = plt_S,
         I         = plt_I,
         R         = plt_R,
         incidence = plt_inc,
         stop("Unsupported 'which'."))
}
