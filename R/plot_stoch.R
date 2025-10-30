#' @title Plot single-population stochastic SIRS outputs
#'
#' @description
#' Visualise **single-population** stochastic SIRS results from a simulation
#' object that contains time series **proportions** (`S`, `I`, `R`) and
#' optionally daily **cases**. Supports four views via `which`:
#' - `"SIR"`: S/I/R mean lines with optional uncertainty ribbons
#' - `"incidence"`: daily cases (mean line + optional ribbon)
#' - `"overlay"`: S/I/R (primary y-axis) with scaled incidence on a secondary axis
#' - `"both"`: two-panel layout with S/I/R (left) and incidence (right)
#'
#' @param sim A single-population stochastic simulation list; see Details for
#'   required components and shapes.
#' @param which Character; one of `c("both","overlay","SIR","incidence")`.
#'   Controls the panel/layout to return.
#' @param probs Numeric; either length-2 `(lower, upper)` quantiles or a single
#'   central mass in `(0,1)` (e.g., `0.95`). Controls uncertainty ribbons.
#' @param per_million Logical; if `TRUE`, scale incidence to per-million using
#'   `sim$params$pop` or `sim$pop` when available. Default `FALSE`.
#' @param sir_cols Named colours for S/I/R (case-insensitive names accepted).
#' @param inc_col Colour for incidence line/ribbon. Default `"grey40"`.
#' @param base_size Numeric base font size for the ggplot theme. Default `12`.
#' @param show_bands Logical; show quantile ribbons when available. Default `TRUE`.
#'
#' @return
#' - For `which = "SIR"`, `"incidence"`, or `"overlay"`: a **ggplot** object.
#' - For `which = "both"`: a **patchwork** object if {patchwork} is installed,
#'   otherwise a list with elements `left` and `right` (two ggplots).
#'
#' @examples
#' \dontrun{
#' # Assume `sim` is from a single-pop stochastic simulator, e.g. simulate_sirs_stoch()
#' p1 <- plot_stoch(sim, which = "SIR", probs = 0.9)
#' p2 <- plot_stoch(sim, which = "incidence", per_million = TRUE)
#' p3 <- plot_stoch(sim, which = "overlay", show_bands = TRUE)
#' p4 <- plot_stoch(sim, which = "both")
#' }
#'
#' @export

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
