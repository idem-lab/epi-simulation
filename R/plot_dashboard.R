# ============================================================
# dashboard.R — ggplot list + arranger (no contact matrix)
# Fixed SIR colours; det/stoch-aware; P=1 friendly; prints parameters
# Auto-coerces sim shapes; returns a named list of ggplots
# ============================================================
# Overview
# - This file builds a *dashboard* of ggplots for SIRS simulations and returns
#   them as a named list (so you can arrange them with patchwork).
# - It auto-detects deterministic (matrices) vs stochastic (arrays) outputs
#   and gracefully handles single- and multi-population cases.
# - It also prints a compact parameter summary to the console and includes a
#   "Parameters" plot panel for static reporting.
#
# Expected shapes
# - Deterministic: S/I/R/Incidence as [T x P] matrices (P can be 1)
# - Stochastic:   proportions [T x sims x P x state], cases [T x sims x P]
# - Time vector:  sim$time of length T
#
# Key panels returned
#   S, I, R          — ribbon (stoch) + mean +/or deterministic line(s)
#   incidence        — incidence per-pop; optional per-million scaling
#   beta             — β(t) as vector or matrix
#   peak_I           — (TODO) peak infected per pop (det) or mean + quantile bars (stoch)
#   stoch_mag        — (TODO) average ribbon width for I (how wide the uncertainty is)
#   det_vs_stoch     — overlay: deterministic vs stochastic (I or incidence)
#   sir_basic        — single-pop S/I/R in one axes when P = 1
#   params           — static text panel of parameters
#
# Usage (minimal)
#   ribbon_probs <- c(0.025, 0.975)
#   plots <- plot_dashboard(sim, probs = ribbon_probs, per_million = TRUE)
#   arrange_dashboard(plots[c("params","I","incidence","beta","det_vs_stoch")], c(2,3))

#' @importFrom stats quantile
NULL

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------- helpers & palettes ----------
.pop_cols <- function(P) {
  c("#3B528B", "#5DC863", "#FDE725", "#21918C", "#440154", "#F8961E", "#58A4B0")[seq_len(min(P,7))]
}
STATE_COLS <- c(S = "#0072B2", I = "#CC79A7", R = "#009E73")

.alpha <- function(cols, a = 0.25) {
  if (length(cols)==1L) grDevices::adjustcolor(cols, a) else vapply(cols, grDevices::adjustcolor, "", alpha.f=a)
}
.not_blank_df <- function(df) !is.null(df) && is.data.frame(df) && nrow(df)>0

.normalize_pop_names <- function(P, pop_names=NULL){
  if (is.null(pop_names)) return(paste0("Pop ", seq_len(P)))
  if (length(pop_names)==P) return(pop_names)
  warning(sprintf("Length of pop_names (%d) != P (%d); using defaults.", length(pop_names), P))
  paste0("Pop ", seq_len(P))
}

# Normalize probs; default to 95% as a vector c(0.025, 0.975).
# (Still accepts a single number p and converts to c(1-p, p) if passed.)
.norm_probs <- function(probs){
  if (missing(probs) || is.null(probs)) return(c(0.025, 0.975))
  stopifnot(is.numeric(probs), all(is.finite(probs)))
  if (length(probs)==1L) {
    p <- as.numeric(probs); stopifnot(p>0, p<1)
    return(c(1-p, p))
  }
  stopifnot(length(probs)>=2)
  sort(probs[1:2])
}

# ---------- AUTO-COERCION ----------
.to_mat_lenT <- function(x,T){
  if (is.null(x)) return(NULL)
  if (is.matrix(x)){ stopifnot(nrow(x)==T); return(x) }
  x <- as.numeric(x); stopifnot(length(x)==T); matrix(x, ncol=1)
}

.ensure_4d_proportions <- function(prop,T){
  if (is.null(prop)) return(NULL)
  d <- dim(prop)
  if (length(d)==4L){
    stopifnot(d[1]==T)
    dn <- dimnames(prop)
    if (is.null(dn) || is.null(dn[[4]])){
      k <- d[4]; labs <- c("S","I","R","X","Y","Z")[seq_len(k)]
      dimnames(prop) <- if (is.null(dn)) list(NULL,NULL,NULL,labs) else { dn[[4]] <- labs; dn }
    }
    return(prop)
  }
  if (length(d)==3L){
    stopifnot(d[1]==T)
    k <- d[3]; labs <- dimnames(prop); st <- if (!is.null(labs)) labs[[3]] else NULL
    if (is.null(st)) st <- c("S","I","R")[seq_len(k)]
    array(prop, dim=c(d[1], d[2], 1L, k), dimnames=list(NULL,NULL,NULL,st))
  } else stop("`proportions` must be 3D [T x sims x state] or 4D [T x sims x pop x state].")
}

.ensure_3d_cases <- function(cases,T){
  if (is.null(cases)) return(NULL)
  d <- dim(cases)
  if (length(d)==3L){ stopifnot(d[1]==T); return(cases) }
  if (length(d)==2L){ stopifnot(d[1]==T); return(array(cases, dim=c(d[1], d[2], 1L))) }
  stop("`cases` must be 2D [T x sims] or 3D [T x sims x pop].")
}

.coerce_sim_for_dashboard <- function(sim){
  stopifnot(!is.null(sim$time))
  T <- length(sim$time)
  S <- .to_mat_lenT(sim$S,T); I <- .to_mat_lenT(sim$I,T); R <- .to_mat_lenT(sim$R,T); Inc <- .to_mat_lenT(sim$incidence,T)
  prop <- .ensure_4d_proportions(sim$proportions,T)
  cases<- .ensure_3d_cases(sim$cases,T)
  P_det <- suppressWarnings(max(c(ncol(S), ncol(I), ncol(R), ncol(Inc)), na.rm=TRUE)); if (!is.finite(P_det)) P_det <- NA_integer_
  P_sto <- if (!is.null(prop)) dim(prop)[3] else NA_integer_
  if (!is.na(P_sto) && !is.na(P_det) && P_det!=P_sto){
    if (P_det==1 && P_sto>1){
      if (!is.null(S))   S   <- matrix(S[,1],   T, P_sto)
      if (!is.null(I))   I   <- matrix(I[,1],   T, P_sto)
      if (!is.null(R))   R   <- matrix(R[,1],   T, P_sto)
      if (!is.null(Inc)) Inc <- matrix(Inc[,1], T, P_sto)
      P_det <- P_sto
    } else if (P_sto==1 && P_det>1){
      prop  <- array(prop,  dim=c(T, dim(prop)[2],  P_det, dim(prop)[4])); for (p in seq_len(P_det)) prop[,,p,] <- prop[,,1,]
      cases <- array(cases, dim=c(T, dim(cases)[2], P_det));              for (p in seq_len(P_det)) cases[,,p]   <- cases[,,1]
      P_sto <- P_det
    } else stop(sprintf("Population mismatch (det P=%s vs stoch P=%s).", P_det, P_sto))
  }
  P <- if (!is.na(P_det)) P_det else if (!is.na(P_sto)) P_sto else 1L
  params <- sim$params %||% list()
  params$pop_vec   <- if (is.null(params$pop_vec) || length(params$pop_vec)!=P) rep(1,P) else params$pop_vec
  params$pop_names <- .normalize_pop_names(P, params$pop_names)
  params$beta      <- params$beta %||% NULL
  list(time=sim$time, S=S, I=I, R=R, incidence=Inc, proportions=prop, cases=cases, params=params)
}

# ---------- LONG-FORM melts ----------
.melt_mat <- function(mat,time,value,pop_names=NULL){
  if (is.null(mat)) return(NULL)
  stopifnot(is.matrix(mat), nrow(mat)==length(time))
  T <- nrow(mat); P <- ncol(mat)
  tibble::tibble(
    time = rep(time, times=P),
    pop  = factor(rep(seq_len(P), each=T), labels=pop_names %||% paste0("Pop ", seq_len(P))),
    !!value := as.vector(mat)
  )
}

# Default 95% band as vector c(0.025, 0.975)
.melt_band <- function(arr,time,probs=c(0.025, 0.975),per_million=FALSE,pop_size=NULL,pop_names=NULL){
  probs <- .norm_probs(probs)
  stopifnot(length(dim(arr))==3, dim(arr)[1]==length(time))
  T <- dim(arr)[1]; P <- dim(arr)[3]
  qL <- apply(arr, c(1,3), stats::quantile, probs[1], na.rm=TRUE)
  qU <- apply(arr, c(1,3), stats::quantile, probs[2], na.rm=TRUE)
  mu <- apply(arr, c(1,3), mean, na.rm=TRUE)
  if (per_million && !is.null(pop_size)) {
    scale <- matrix(pop_size/1e6, nrow=T, ncol=P, byrow=TRUE)
    mu <- mu/scale; qL <- qL/scale; qU <- qU/scale
  }
  tibble::tibble(
    time = rep(time, times=P),
    pop  = factor(rep(seq_len(P), each=T), labels=pop_names %||% paste0("Pop ", seq_len(P))),
    mean = as.vector(mu),
    qL   = as.vector(qL),
    qU   = as.vector(qU)
  )
}

# ---- console helper ----
.print_params <- function(lines_vec){
  if (requireNamespace("cli", quietly=TRUE)) { cli::cli_h2("Simulation parameters"); cli::cli_ul(lines_vec) }
  else { cat("\n=== Simulation parameters ===\n"); cat(paste0(" - ", lines_vec, "\n"), sep="") }
}

# ---------- panel builders ----------
.build_SIR_panel <- function(state_name, det_df, band_df, ylab, P, facet_cols=2,
                             group_style=c("facet","combined"), cols=NULL, show_bands=TRUE){
  stopifnot(state_name %in% c("S","I","R"))
  group_style <- match.arg(group_style)
  state_col <- unname(STATE_COLS[[state_name]])
  single_pop <- (P==1)
  
  if (group_style=="combined" && !single_pop){
    p <- ggplot2::ggplot()
    if (.not_blank_df(det_df)) {
      p <- p + ggplot2::geom_line(data=det_df, ggplot2::aes(time, .data[[state_name]], color=pop), linewidth=0.9) +
        ggplot2::scale_color_manual(values=cols, drop=FALSE, name="Population")
    } else if (.not_blank_df(band_df)) {
      p <- p + ggplot2::geom_line(data=band_df, ggplot2::aes(time, mean, color=pop), linewidth=0.9) +
        ggplot2::scale_color_manual(values=cols, drop=FALSE, name="Population")
    }
    return(p + ggplot2::labs(title=paste0(state_name,"(t)"), x="Day", y=ylab) +
             ggplot2::theme_minimal(base_size=11) +
             ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)))
  }
  
  p <- ggplot2::ggplot()
  if (.not_blank_df(band_df)) {
    if (show_bands) {
      p <- p + ggplot2::geom_ribbon(data=band_df, ggplot2::aes(time, ymin=qL, ymax=qU),
                                    inherit.aes=FALSE, fill=grDevices::adjustcolor(state_col, 0.25))
    }
    p <- p + ggplot2::geom_line(data=band_df, ggplot2::aes(time, mean, group=pop),
                                inherit.aes=FALSE, linewidth=0.9, color=state_col)
  }
  if (.not_blank_df(det_df)) {
    p <- p + ggplot2::geom_line(data=det_df, ggplot2::aes(time, .data[[state_name]]),
                                inherit.aes=FALSE, linewidth=0.9, color=state_col)
  }
  if (!single_pop && (.not_blank_df(band_df) || .not_blank_df(det_df))) {
    p <- p + ggplot2::facet_wrap(~pop, ncol=facet_cols)
  }
  p + ggplot2::labs(title=paste0(state_name,"(t)"), x="Day", y=ylab) +
    ggplot2::theme_minimal(base_size=11) +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}

.build_overlay_detstoch <- function(state=c("I","incidence"), det_df, band_df, P, cols,
                                    det_col="#0072B2", title="Det vs Stoch", facet_cols=2, show_bands=TRUE){
  state <- match.arg(state)
  single_pop <- (P==1)
  yvar <- if (state=="I") "I" else "cases"
  if ((is.null(band_df) || nrow(band_df)==0) && (is.null(det_df) || nrow(det_df)==0)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title=paste(title,"(no data)")))
  }
  if (single_pop){
    p <- ggplot2::ggplot()
    if (!is.null(band_df) && nrow(band_df)>0){
      if (show_bands){
        p <- p + ggplot2::geom_ribbon(data=band_df, ggplot2::aes(time, ymin=qL, ymax=qU, fill="Stochastic band"), alpha=0.25)
      }
      p <- p + ggplot2::geom_line(data=band_df, ggplot2::aes(time, mean, color="Stochastic mean"), linewidth=0.9)
    }
    if (!is.null(det_df) && nrow(det_df)>0){
      p <- p + ggplot2::geom_line(data=det_df, ggplot2::aes(time, .data[[yvar]], color="Deterministic"), linewidth=0.9)
    }
    return(p +
             ggplot2::scale_fill_manual(values=c("Stochastic band"=grDevices::adjustcolor("grey50",0.25)), name=NULL) +
             ggplot2::scale_color_manual(values=c("Stochastic mean"="black","Deterministic"=det_col), name=NULL) +
             ggplot2::labs(title=title, x="Day", y=if (state=="I") "I proportion" else "Daily cases") +
             ggplot2::theme_minimal(base_size=11) +
             ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)))
  } else {
    p <- ggplot2::ggplot()
    if (!is.null(band_df) && nrow(band_df)>0){
      if (show_bands){
        p <- p + ggplot2::geom_ribbon(data=band_df, ggplot2::aes(time, ymin=qL, ymax=qU),
                                      fill=grDevices::adjustcolor("grey50",0.25))
      }
      p <- p + ggplot2::geom_line(data=band_df, ggplot2::aes(time, mean), color="black", linewidth=0.7)
    }
    if (!is.null(det_df) && nrow(det_df)>0){
      p <- p + ggplot2::geom_line(data=det_df, ggplot2::aes(time, .data[[yvar]], color=pop), linewidth=0.9)
    }
    if (.not_blank_df(band_df) || .not_blank_df(det_df)) p <- p + ggplot2::facet_wrap(~pop, ncol=facet_cols)
    p + ggplot2::scale_color_manual(values=cols, drop=FALSE, name="Population") +
      ggplot2::labs(title=title, x="Day", y=if (state=="I") "I proportion" else "Daily cases") +
      ggplot2::theme_minimal(base_size=11) +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
  }
}

# ---------- core API ----------
plot_dashboard <- function(sim,
                           probs = c(0.025, 0.975),  # default 95% as vector
                           per_million = FALSE,
                           title = "Dashboard",
                           overlay_state = c("I","incidence"),
                           overlay_probs = NULL,
                           overlay_det_col = "#0072B2",
                           overlay_title   = "Det vs Stoch (overlay)",
                           group_style = c("facet","combined"),
                           pop_cols = NULL,
                           show_bands = TRUE){
  if (!requireNamespace("ggplot2", quietly=TRUE)) stop("plot_dashboard() needs {ggplot2}.")
  group_style  <- match.arg(group_style)
  overlay_state<- match.arg(overlay_state)
  
  # normalize probs (accept c(L,U) *or* single number p)
  probs <- .norm_probs(probs)
  
  sim  <- .coerce_sim_for_dashboard(sim)
  time <- sim$time
  has_det <- !is.null(sim$I) && is.matrix(sim$I)
  has_sto <- !is.null(sim$proportions) && length(dim(sim$proportions))==4
  if (!has_det && !has_sto) stop("Unsupported `sim` shape: need det matrices or stoch arrays.")
  
  P         <- if (has_det) ncol(sim$I) else dim(sim$proportions)[3]
  pop_names <- .normalize_pop_names(P, sim$params$pop_names)
  pops      <- sim$params$pop_vec %||% rep(1,P)
  cols <- if (!is.null(pop_cols)) { stopifnot(length(pop_cols)>=P); setNames(unname(pop_cols)[seq_len(P)], pop_names) }
  else { setNames(.pop_cols(P), pop_names) }
  
  inc_lab <- if (per_million) "Daily cases (per million)" else "Daily cases"
  s_lab <- "S proportion"; i_lab <- "I proportion"; r_lab <- "R proportion"
  
  S_det <- if (has_det) .melt_mat(sim$S, time, "S", pop_names) else NULL
  I_det <- if (has_det) .melt_mat(sim$I, time, "I", pop_names) else NULL
  R_det <- if (has_det) .melt_mat(sim$R, time, "R", pop_names) else NULL
  Inc_det <- NULL
  if (has_det) {
    inc_mat <- sim$incidence; if (per_million) inc_mat <- sweep(inc_mat,2,pops/1e6,"/")
    Inc_det <- .melt_mat(inc_mat, time, "cases", pop_names)
  }
  
  I_band <- Inc_band <- S_band <- R_band <- NULL
  if (has_sto) {
    .get_state_3d <- function(prop4, state){
      d <- dim(prop4); dn4 <- dimnames(prop4); stn <- if (is.null(dn4)) NULL else dn4[[4]]
      idx <- if (!is.null(stn) && state %in% stn) which(stn==state)[1] else { i <- switch(state,S=1L,I=2L,R=3L,1L); if (i>d[4]) 1L else i }
      array(prop4[,,,idx, drop=FALSE], dim=d[1:3])
    }
    S_arr <- .get_state_3d(sim$proportions,"S")
    I_arr <- .get_state_3d(sim$proportions,"I")
    R_arr <- .get_state_3d(sim$proportions,"R")
    cases_a <- sim$cases; if (length(dim(cases_a))==2L) cases_a <- array(cases_a, dim=c(dim(cases_a),1L))
    
    S_band   <- .melt_band(S_arr,   time, probs, FALSE, NULL, pop_names)
    I_band   <- .melt_band(I_arr,   time, probs, FALSE, NULL, pop_names)
    R_band   <- .melt_band(R_arr,   time, probs, FALSE, NULL, pop_names)
    Inc_band <- .melt_band(cases_a, time, probs, per_million, pops, pop_names)
  }
  
  # ---- parameters panel (static text) ----
  p <- sim$params
  params_plot <- {
    line <- function(lbl,val) sprintf("%s: %s", lbl, val)
    beta_str <- if (is.null(p$beta)) "—" else if (is.matrix(p$beta)) "[matrix]" else if (length(p$beta)==length(time)) "[vector]" else paste(p$beta, collapse=",")
    lines_vec <- c(
      line("n_times", length(time)),
      line("P", P),
      line("pop_vec", if (is.null(p$pop_vec)) "—" else paste(p$pop_vec, collapse=", ")),
      line("I_init",  if (is.null(p$I_init))  "—" else paste(p$I_init, collapse=", ")),
      line("beta", beta_str),
      line("gamma", if (is.null(p$gamma)) "—" else signif(p$gamma,4)),
      line("1/gamma (days)", if (is.null(p$gamma)) "—" else sprintf("%.2f", 1/p$gamma)),
      line("omega", if (is.null(p$omega)) "—" else signif(p$omega,4)),
      line("1/omega (days)", if (is.null(p$omega)) "—" else sprintf("%.2f", 1/p$omega)),
      line("epsilon", if (is.null(p$epsilon)) "—" else signif(p$epsilon,4)),
      line("alpha",   if (is.null(p$alpha))   "—" else signif(p$alpha,4)),
      line("n_sims",  if (!is.null(p$n_sims)) p$n_sims else (if (has_det) 1 else "—")),
      line("stochastic", if (!is.null(p$stochastic)) p$stochastic else has_sto)
    )
    .print_params(lines_vec)
    y_start <- 0.92; y_step <- 0.095; yy <- y_start - y_step*seq_along(lines_vec) + y_step
    df <- data.frame(x=0.5, y=yy, label=lines_vec)
    ggplot2::ggplot(df, ggplot2::aes(x, y, label=label)) +
      ggplot2::geom_text(size=3.5, lineheight=1.0, hjust=0.5, vjust=1) +
      ggplot2::annotate("rect", xmin=0, xmax=1, ymin=0, ymax=1, fill=NA, colour="grey60") +
      ggplot2::annotate("text", x=0.5, y=0.98, label="Parameters", fontface="bold", size=4, vjust=1) +
      ggplot2::scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
      ggplot2::scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
      ggplot2::theme_void() + ggplot2::theme(plot.margin=grid::unit(c(6,6,6,6), "pt"))
  }
  
  # basic SIR single-pop
  plt_sir_basic <- NULL
  if (P==1){
    get_SIR <- function(){
      if (has_det && !is.null(sim$S) && !is.null(sim$I) && !is.null(sim$R))
        list(S=sim$S[,1], I=sim$I[,1], R=sim$R[,1])
      else if (has_sto && .not_blank_df(S_band) && .not_blank_df(I_band) && .not_blank_df(R_band)){
        pop1 <- levels(S_band$pop)[1]
        list(S=S_band$mean[S_band$pop==pop1], I=I_band$mean[I_band$pop==pop1], R=R_band$mean[R_band$pop==pop1])
      } else NULL
    }
    sir <- get_SIR()
    if (!is.null(sir)){
      Tlen <- length(time)
      df_basic <- tibble::tibble(
        time  = rep(time, times=3L),
        state = factor(rep(c("S","I","R"), each=Tlen), levels=c("S","I","R")),
        value = c(as.numeric(sir$S), as.numeric(sir$I), as.numeric(sir$R))
      )
      plt_sir_basic <- ggplot2::ggplot(df_basic, ggplot2::aes(time, value, colour=state)) +
        ggplot2::geom_line(linewidth=0.9) +
        ggplot2::scale_colour_manual(values=STATE_COLS, name=NULL) +
        ggplot2::labs(title="S, I, R (single population)", x="Day", y="Proportion") +
        ggplot2::theme_minimal(base_size=11) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
    }
  }
  
  # panels
  plt_S <- .build_SIR_panel("S", S_det, S_band, s_lab, P, group_style=group_style, cols=cols, show_bands=show_bands)
  plt_I <- .build_SIR_panel("I", I_det, I_band, i_lab, P, group_style=group_style, cols=cols, show_bands=show_bands)
  plt_R <- .build_SIR_panel("R", R_det, R_band, r_lab, P, group_style=group_style, cols=cols, show_bands=show_bands)
  
  # incidence
  plt_inc <- {
    single_pop <- (P==1)
    pplt <- ggplot2::ggplot()
    if (.not_blank_df(Inc_band)) {
      if (show_bands) {
        pplt <- pplt + ggplot2::geom_ribbon(ggplot2::aes(time, ymin=qL, ymax=qU, fill=pop), data=Inc_band, alpha=0.25, show.legend=!single_pop)
      }
      pplt <- pplt + ggplot2::geom_line(ggplot2::aes(time, mean, color=pop), data=Inc_band, linewidth=0.9, show.legend=!single_pop)
    }
    if (.not_blank_df(Inc_det)) {
      pplt <- pplt + ggplot2::geom_line(ggplot2::aes(time, cases, color=pop), data=Inc_det, linewidth=0.9, show.legend=!single_pop)
    }
    pplt +
      ggplot2::scale_color_manual(values=cols, drop=FALSE) +
      ggplot2::scale_fill_manual(values=setNames(.alpha(cols,0.25), names(cols)), drop=FALSE) +
      ggplot2::labs(title="Incidence", x="Day", y=inc_lab, color="Population", fill="Population") +
      ggplot2::theme_minimal(base_size=11) +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),
                     legend.position = if (single_pop) "none" else "right")
  }
  
  # beta(t)
  beta_obj <- sim$params$beta
  plt_beta <- {
    if (is.null(beta_obj)) {
      ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title=expression(beta*"(t) unavailable"))
    } else if (is.matrix(beta_obj) && nrow(beta_obj)==length(time) && ncol(beta_obj)==P) {
      beta_long <- .melt_mat(beta_obj, time, "beta", pop_names)
      ggplot2::ggplot(beta_long, ggplot2::aes(time, beta, color=pop)) +
        ggplot2::geom_line(linewidth=0.9) +
        ggplot2::scale_color_manual(values=cols) +
        ggplot2::labs(title=expression(beta*"(t) by population"), x="Day", y=expression(beta), color="Population") +
        ggplot2::theme_minimal(base_size=11) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
    } else if (is.vector(beta_obj) && length(beta_obj)==length(time)) {
      ggplot2::ggplot(tibble::tibble(time=time, beta=as.numeric(beta_obj)), ggplot2::aes(time, beta)) +
        ggplot2::geom_line(linewidth=0.9) +
        ggplot2::labs(title=expression(beta*"(t)"), x="Day", y=expression(beta)) +
        ggplot2::theme_minimal(base_size=11) + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
    } else ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title=expression(beta*"(t) shape mismatch"))
  }
  
  # Det vs Stoch overlay
  plt_det_vs_stoch <- {
    ost <- overlay_state
    q_overlay <- if (is.null(overlay_probs)) probs else .norm_probs(overlay_probs)
    band_df <- if (ost=="I") I_band else Inc_band
    det_df  <- if (ost=="I") I_det  else Inc_det
    if (!is.null(overlay_probs)) {
      if (ost=="I" && !is.null(sim$proportions)) {
        arr <- sim$proportions[,,, "I", drop=FALSE]; arr <- array(arr, dim=dim(arr)[1:3])
        band_df <- .melt_band(arr, time, q_overlay, FALSE, NULL, pop_names)
      } else if (ost=="incidence" && !is.null(sim$cases)) {
        arr <- sim$cases; if (length(dim(arr))==2L) arr <- array(arr, dim=c(dim(arr),1L))
        band_df <- .melt_band(arr, time, q_overlay, per_million, pops, pop_names)
      }
    }
    .build_overlay_detstoch(ost, det_df, band_df, P, cols, overlay_det_col, overlay_title, show_bands=show_bands)
  }
  
  list(
    S = plt_S,
    I = plt_I,
    R = plt_R,
    incidence = plt_inc,
    beta = plt_beta,
    det_vs_stoch = plt_det_vs_stoch,
    SIR = plt_sir_basic,
    params = params_plot
  )
}

# ---------- arranger ----------
arrange_dashboard <- function(plots, layout=c(2,3), collect_legend=TRUE){
  if (!requireNamespace("patchwork", quietly=TRUE)) stop("arrange_dashboard() needs {patchwork}.")
  if (!requireNamespace("ggplot2", quietly=TRUE))  stop("arrange_dashboard() needs {ggplot2}.")
  stopifnot(is.list(plots), length(layout)==2)
  pw <- patchwork::wrap_plots(plots, nrow=layout[1], ncol=layout[2], guides=if (collect_legend) "collect" else "auto", byrow=TRUE)
  if (collect_legend) pw <- pw & ggplot2::theme(legend.position="right")
  pw
}
