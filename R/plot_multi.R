# ===========================
# plot_stoch — single-pop stoch
# ===========================
plot_stoch <- function(
    sim,
    which        = c("both","overlay","SIR","incidence"),
    probs        = 0.95,
    per_million  = FALSE,
    sir_cols     = c(S="#0072B2", I="#CC79A7", R="#009E73"),
    inc_col      = "grey40",
    base_size    = 12,
    show_bands   = TRUE
){
  which <- match.arg(which)
  stopifnot(!is.null(sim$time), !is.null(sim$proportions))
  time <- sim$time
  prop <- sim$proportions
  
  # --- normalize sir_cols to named S/I/R (case-insensitive, handles unnamed/partial) ---
  .norm_sir_cols <- function(cols) {
    nm <- names(cols)
    if (is.null(nm) || any(nm == "")) {
      cols <- as.vector(cols)
      stopifnot(length(cols) >= 3)
      out <- c(S = cols[1], I = cols[2], R = cols[3])
      return(out)
    }
    key <- toupper(nm)
    map <- setNames(cols, key)
    out <- c(
      S = map[["S"]] %||% cols[[1]],
      I = map[["I"]] %||% cols[[2]],
      R = map[["R"]] %||% cols[[3]]
    )
    names(out) <- c("S","I","R")
    out
  }
  `%||%` <- function(x, y) if (is.null(x) || length(x)==0 || is.na(x)) y else x
  sir_cols <- .norm_sir_cols(sir_cols)
  
  # ---- normalize shapes to [T x sims] per state ----
  slice_state <- function(state){
    if (length(dim(prop)) == 4L) {
      stn <- dimnames(prop)[[4]]
      idx <- if (!is.null(stn)) match(state, stn) else switch(state, S=1L, I=2L, R=3L)
      if (is.na(idx)) stop("State '", state, "' not found in dimnames(proportions).")
      a <- prop[,,, idx, drop = FALSE]; a[,,1,1, drop = TRUE]
    } else if (length(dim(prop)) == 3L) {
      stn <- dimnames(prop)[[3]]
      idx <- if (!is.null(stn)) match(state, stn) else switch(state, S=1L, I=2L, R=3L)
      if (is.na(idx)) stop("State '", state, "' not found in dimnames(proportions).")
      a <- prop[,, idx, drop = FALSE]; a[,,1, drop = TRUE]
    } else stop("proportions must be 3D or 4D.")
  }
  
  cases <- sim$cases
  if (!is.null(cases)) {
    if (length(dim(cases)) == 3L) cases <- cases[,,1, drop = TRUE]  # [T x sims]
    if (!is.null(dim(cases)) && length(dim(cases)) != 2L)
      stop("cases must be [T x sims] (or [T x sims x 1]).")
  }
  
  # ---- quantiles + means per time ----
  norm_probs <- function(p){
    stopifnot(is.numeric(p), all(is.finite(p)))
    if (length(p) == 1L) { p <- as.numeric(p); stopifnot(p>0 && p<1); return(c(1-p, p)) }
    p[1:2]
  }
  q <- norm_probs(probs)
  
  q_df <- function(A, name){
    m  <- rowMeans(A, na.rm = TRUE)
    qL <- apply(A, 1L, stats::quantile, probs = q[1], na.rm = TRUE)
    qU <- apply(A, 1L, stats::quantile, probs = q[2], na.rm = TRUE)
    data.frame(time = time, mean = m, qL = qL, qU = qU, state = name, check.names = FALSE)
  }
  
  df_S <- q_df(slice_state("S"), "S")
  df_I <- q_df(slice_state("I"), "I")
  df_R <- q_df(slice_state("R"), "R")
  
  pop <- if (!is.null(sim$params$pop)) sim$params$pop else sim$pop
  scale_cases <- function(x) {
    if (per_million && is.finite(pop) && pop > 0) return(1e6 * x / pop)
    x
  }
  inc_label <- if (per_million && is.finite(pop) && pop > 0) "Daily cases (per million)" else "Daily cases"
  
  df_inc <- NULL
  if (!is.null(cases)) {
    A <- scale_cases(cases)
    m  <- rowMeans(A, na.rm = TRUE)
    qL <- apply(A, 1L, stats::quantile, probs = q[1], na.rm = TRUE)
    qU <- apply(A, 1L, stats::quantile, probs = q[2], na.rm = TRUE)
    df_inc <- data.frame(time = time, mean = m, qL = qL, qU = qU)
  }
  
  suppressPackageStartupMessages({ requireNamespace("ggplot2") })
  library(ggplot2)
  
  # ---------- builders ----------
  build_sir <- function(){
    df_all <- rbind(df_S, df_I, df_R)
    df_all$state <- factor(df_all$state, levels = c("S","I","R"))
    col_map <- c(S=sir_cols[["S"]], I=sir_cols[["I"]], R=sir_cols[["R"]])
    
    p <- ggplot(df_all, aes(time, mean, color = state, fill = state))
    if (show_bands) p <- p + geom_ribbon(aes(ymin = qL, ymax = qU), alpha = 0.25, colour = NA)
    p + geom_line(linewidth = 1) +
      scale_color_manual(values = col_map, breaks = c("S","I","R"), name = NULL) +
      scale_fill_manual(values  = col_map, breaks = c("S","I","R"), guide = "none") +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Day", y = "Proportion", title = "S, I, R (stochastic)") +
      theme_minimal(base_size = base_size) +
      theme(panel.grid.minor = element_blank(), legend.position = "top")
  }
  
  build_inc <- function(){
    if (is.null(df_inc)) {
      return(ggplot() + theme_void() + ggtitle("Incidence (cases missing)"))
    }
    p <- ggplot(df_inc, aes(time, mean))
    if (show_bands) p <- p + geom_ribbon(aes(ymin = qL, ymax = qU), fill = scales::alpha(inc_col, 0.25))
    p + geom_line(color = inc_col, linewidth = 1) +
      labs(x = "Day", y = inc_label, title = "Daily incidence (stochastic)") +
      theme_minimal(base_size = base_size) +
      theme(panel.grid.minor = element_blank())
  }
  
  build_overlay <- function(){
    p <- build_sir()
    if (!is.null(df_inc)) {
      inc_max <- max(df_inc$qU, na.rm = TRUE)
      k <- if (is.finite(inc_max) && inc_max > 0) 1 / inc_max else 1
      p <- p +
        { if (show_bands) geom_ribbon(data = transform(df_inc, mean = mean * k, qL = qL * k, qU = qU * k),
                                      aes(time, ymin = qL, ymax = qU),
                                      inherit.aes = FALSE,
                                      fill = scales::alpha(inc_col, 0.25)) } +
        geom_line(data = transform(df_inc, mean = mean * k),
                  aes(time, y = mean),
                  inherit.aes = FALSE,
                  color = inc_col, linewidth = 1) +
        scale_y_continuous(
          name = "Proportion",
          limits = c(0, 1),
          sec.axis = sec_axis(~ . / k, name = inc_label)
        )
    } else {
      p <- p + ggtitle("S, I, R (stochastic) — incidence not available")
    }
    p + labs(title = "S, I, R + incidence (stochastic)")
  }
  
  # ---------- dispatch ----------
  if (which == "SIR")       return(build_sir())
  if (which == "incidence") return(build_inc())
  if (which == "overlay")   return(build_overlay())
  
  # which == "both"
  p_left  <- build_sir()
  p_right <- build_inc()
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return((p_left + p_right) + patchwork::plot_layout(ncol = 2, widths = c(2, 1)))
  } else {
    warning("Package {patchwork} not installed; returning a list of two plots.")
    return(list(left = p_left, right = p_right))
  }
}


# ============================================================
# plot_multi — MULTI-POP (det + stoch in one)
# S, I, R: facet = state colours | combined = pop colours
# Incidence: always combined (pop colours)
# ============================================================

plot_multi <- function(
    sim,
    probs       = c(0.025, 0.975),   # default 95% as vector or single p
    per_million = FALSE,
    which       = c("S", "I", "R", "incidence"),
    group_style = c("facet", "combined"),
    pop_cols    = NULL,
    show_bands  = TRUE
){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_multi() needs {ggplot2}.")
  }
  
  which       <- match.arg(which)
  group_style <- match.arg(group_style)
  
  # normalize probs (accept c(L,U) *or* single number p)
  probs <- .norm_probs(probs)
  
  sim  <- .coerce_sim_for_dashboard(sim)
  time <- sim$time
  
  has_det <- !is.null(sim$I) && is.matrix(sim$I)
  has_sto <- !is.null(sim$proportions) && length(dim(sim$proportions)) == 4L
  if (!has_det && !has_sto) {
    stop("Unsupported `sim` shape: need deterministic matrices or stochastic arrays.")
  }
  
  # ---- basic shape checks to avoid 'S/I/R must align with time' ----
  if (has_det) {
    Tlen <- length(time)
    for (nm in c("S","I","R","incidence")) {
      if (!is.null(sim[[nm]])) {
        if (!is.matrix(sim[[nm]])) stop(nm, " must be a [T x P] matrix.")
        if (nrow(sim[[nm]]) != Tlen) stop("Time length and ", nm, " rows differ (", Tlen, " vs ", nrow(sim[[nm]]), ").")
      }
    }
  } else {
    # stochastic: proportions [T x sims x P x state], cases [T x sims x P]
    Tlen <- length(time)
    if (dim(sim$proportions)[1] != Tlen) {
      stop("Time length and proportions first dim differ (", Tlen, " vs ", dim(sim$proportions)[1], ").")
    }
    if (!is.null(sim$cases) && dim(sim$cases)[1] != Tlen) {
      stop("Time length and cases first dim differ (", Tlen, " vs ", dim(sim$cases)[1], ").")
    }
  }
  
  P         <- if (has_det) ncol(sim$I) else dim(sim$proportions)[3]
  pop_names <- .normalize_pop_names(P, sim$params$pop_names)
  pops      <- sim$params$pop_vec %||% rep(1, P)
  
  # ---- palette: prefer dashboard helper; else pop_cols -> my_pop_cols -> fallback ----
  if (exists(".reconcile_pop_cols", mode = "function")) {
    cols <- .reconcile_pop_cols(P, pop_names, pop_cols)
  } else {
    base_cols <- if (is.null(pop_cols)) {
      gcols <- get0("my_pop_cols", ifnotfound = NULL, inherits = TRUE)
      if (is.null(gcols)) .pop_cols(P) else unname(as.character(gcols))
    } else {
      unname(as.character(pop_cols))
    }
    if (length(base_cols) < P) base_cols <- rep(base_cols, length.out = P)
    if (length(base_cols) > P) base_cols <- base_cols[seq_len(P)]
    cols <- stats::setNames(base_cols, pop_names)
  }
  
  inc_lab <- if (per_million) "Daily cases (per million)" else "Daily cases"
  s_lab <- "S proportion"; i_lab <- "I proportion"; r_lab <- "R proportion"
  
  # ---- deterministic melts ----
  S_det <- I_det <- R_det <- Inc_det <- NULL
  if (has_det) {
    S_det <- .melt_mat(sim$S, time, "S", pop_names)
    I_det <- .melt_mat(sim$I, time, "I", pop_names)
    R_det <- .melt_mat(sim$R, time, "R", pop_names)
    inc_mat <- sim$incidence
    if (!is.null(inc_mat)) {
      if (per_million) inc_mat <- sweep(inc_mat, 2, pops / 1e6, "/")
      Inc_det <- .melt_mat(inc_mat, time, "cases", pop_names)
    }
  }
  
  # ---- stochastic bands ----
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
  
  # ---- panel builders ----
  plt_S <- .build_SIR_panel("S", S_det, S_band, s_lab, P,
                            group_style = group_style, cols = cols, show_bands = show_bands)
  plt_I <- .build_SIR_panel("I", I_det, I_band, i_lab, P,
                            group_style = group_style, cols = cols, show_bands = show_bands)
  plt_R <- .build_SIR_panel("R", R_det, R_band, r_lab, P,
                            group_style = group_style, cols = cols, show_bands = show_bands)
  
  # incidence
  plt_inc <- {
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
    pplt +
      ggplot2::scale_color_manual(values = cols, drop = FALSE) +
      ggplot2::scale_fill_manual(values = stats::setNames(.alpha(cols, 0.25), names(cols)), drop = FALSE) +
      ggplot2::labs(title = "Daily Incidence", x = "Day", y = inc_lab, color = "Population", fill = "Population") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = if (single_pop) "none" else "right"
      )
  }
  
  # ---- return the requested panel ----
  switch(which,
         S         = plt_S,
         I         = plt_I,
         R         = plt_R,
         incidence = plt_inc,
         stop("Unsupported 'which'."))
}
