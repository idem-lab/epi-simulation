#' Plot multi-population SIRS outputs (deterministic + stochastic)
#'
#' @description
#' Visualise SIRS simulation results for **multiple populations** from either
#' deterministic runs (`[T x P]` matrices) or stochastic runs
#' (`proportions [T x sims x P x {S,I,R}]`, `cases [T x sims x P]`).
#' Supports:
#' - **S / I / R panels**:  
#'   - `group_style = "facet"` → facet by **state**, use fixed S/I/R colours.  
#'   - `group_style = "combined"` → one panel, colours map to **populations**.  
#' - **Incidence panel** (always combined): mean line per population with optional
#'   quantile ribbon across simulations.
#'
#' @param sim A simulation result list. Deterministic inputs should provide
#'   `[T x P]` matrices for `S`, `I`, `R`, and optionally `incidence`, plus
#'   a `time` vector. Stochastic inputs should provide `proportions`
#'   `[T x sims x P x {S,I,R}]`, optionally `cases [T x sims x P]`, and `time`.
#' @param probs Length-2 numeric vector of lower/upper quantiles for the ribbon,
#'   or a single number in `(0,1)` interpreted as a symmetric central mass
#'   (e.g., `0.95` → `c(0.025, 0.975)`). Default `c(0.025, 0.975)`.
#' @param per_million Logical; if `TRUE`, incidence is scaled to per-million
#'   using population sizes (if available). Default `FALSE`.
#' @param which Character; which panel to return. One of
#'   `c("S","I","R","incidence")`.
#' @param group_style Character; how to group colours for S/I/R panels:
#'   `"facet"` (facet by state; fixed S/I/R colours) or `"combined"`
#'   (single panel; colours map to populations). Ignored for `"incidence"`.
#' @param pop_cols Optional character vector of colours to use for populations.
#'   Recycled or truncated to match the number of populations if needed.
#' @param show_bands Logical; show quantile ribbons for stochastic inputs.
#'   Has no effect for deterministic-only inputs. Default `TRUE`.
#'
#' @return A **ggplot** object corresponding to the selected `which` panel.
#'
#' @examples
#' \dontrun{
#' # Deterministic multi-pop (toy)
#' sim_det <- simulate_sirs_multi_det(
#'   n_times = 120,
#'   pop_vec = c(8e4, 6e4),
#'   I_init  = c(8, 6),
#'   beta_mat = matrix(0.16, nrow = 120, ncol = 2),
#'   gamma = 1/7, omega = 1/60
#' )
#' p1 <- plot_multi(sim_det, which = "I", group_style = "combined")
#' p2 <- plot_multi(sim_det, which = "incidence", per_million = TRUE)
#'
#' # Stochastic multi-pop (toy)
#' sim_sto <- simulate_sirs_multi_stoch(
#'   n_times = 120, pop_vec = c(8e4, 6e4), I_init = c(8, 6),
#'   beta_mat = matrix(0.16, nrow = 120, ncol = 2),
#'   gamma = 1/7, omega = 1/60, n_sims = 50, seed = 42
#' )
#' p3 <- plot_multi(sim_sto, which = "S", group_style = "facet", probs = c(0.1, 0.9))
#' p4 <- plot_multi(sim_sto, which = "incidence", show_bands = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_color_manual
#'   scale_fill_manual labs theme_minimal element_text
#' @importFrom stats setNames
#' @export

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
