#' @title Plot SIR trajectories and/or incidence diagnostics
#'
#' @description
#' Produce diagnostic plots from a single-population SIR/SIRS simulation:
#' - `"sir"`: S / I / R proportions over time
#' - `"incidence"`: daily new infections (counts)
#' - `"both_side"`: S/I/R and incidence side by side
#' - `"overlay"`: S/I/R plus incidence on a secondary y-axis
#'
#' @param sim A list-like simulation object with numeric vectors:
#'   `time`, `S`, `I`, `R`, and `incidence` (all same length).
#' @param which Character. One of
#'   `c("both_side","overlay","sir","incidence")`.
#'
#' @return A ggplot (for "sir", "incidence", "overlay") or a patchwork object
#'   (for "both_side"). This makes it easy to print, save with `ggsave()`,
#'   or combine with other plots. No invisible(NULL).
#'
#' @examples
#' # Minimal reproducible example
#' set.seed(1)
#' t <- 1:120
#' S <- pmax(0, 1 - 0.6*(1 - exp(-t/60)))
#' I <- pmax(0, 0.15*exp(-(t-40)^2/800))
#' R <- pmax(0, 1 - S - I)
#' inc <- pmax(0, round(diff(c(0, I))*1000 + rnorm(length(t), 0, 2)))
#' sim <- list(time = t, S = S, I = I, R = R, incidence = inc)
#'
#' p1 <- plot_sirs(sim, which = "sir")
#' p2 <- plot_sirs(sim, which = "incidence")
#' p3 <- plot_sirs(sim, which = "both_side")
#' p4 <- plot_sirs(sim, which = "overlay")
#'
#' @export
plot_sirs <- function(
    sim,
    which = c("both_side","overlay","sir","incidence")
) {
  # ---- pick mode ----
  which <- match.arg(which)
  
  # ---- basic checks ----
  need <- c("time","S","I","R","incidence")
  miss <- setdiff(need, names(sim))
  if (length(miss)) stop("sim is missing: ", paste(miss, collapse=", "))
  
  # coerce to data.frames we can ggplot
  df_sir <- data.frame(
    time = sim$time,
    S    = sim$S,
    I    = sim$I,
    R    = sim$R
  )
  
  df_inc <- data.frame(
    time      = sim$time,
    incidence = sim$incidence
  )
  
  # long format for S/I/R (for nice legend)
  df_long <- tidyr::pivot_longer(
    df_sir,
    cols = c("S","I","R"),
    names_to = "state",
    values_to = "value"
  )
  
  # fixed colours for S / I / R
  state_cols <- c(
    S = "#0072B2",  # blue
    I = "#CC79A7",  # magenta
    R = "#009E73"   # green
  )
  
  # ---------------------------
  # base plot pieces
  # ---------------------------
  
  p_sir <- ggplot2::ggplot(df_long,
                           ggplot2::aes(x = time, y = value, colour = state)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(values = state_cols, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0,1), name = "proportion") +
    ggplot2::labs(
      x = "day",
      title = "S, I, R"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  p_inc <- ggplot2::ggplot(df_inc,
                           ggplot2::aes(x = time, y = incidence)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "day",
      y = "new infections (count)",
      title = "Daily incidence"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  # ---------------------------
  # overlay plot with 2nd axis
  # ---------------------------
  # We'll rescale incidence to [0,1] range so it can sit on top of S/I/R.
  # Then we add a secondary axis that maps it back to counts.
  max_inc <- max(df_inc$incidence, na.rm = TRUE)
  if (!is.finite(max_inc) || max_inc == 0) {
    max_inc <- 1  # avoid divide-by-zero
  }
  
  df_overlay <- dplyr::left_join(
    df_long,
    df_inc,
    by = "time"
  )
  df_overlay$incidence_scaled <- df_overlay$incidence / max_inc
  
  # We'll plot two geoms:
  # - lines for S/I/R (coloured)
  # - line for incidence_scaled (grey/black)
  p_overlay <- ggplot2::ggplot() +
    # S / I / R
    ggplot2::geom_line(
      data = df_long,
      ggplot2::aes(x = time, y = value, colour = state),
      linewidth = 1
    ) +
    # incidence on rescaled axis
    ggplot2::geom_line(
      data = df_inc,
      ggplot2::aes(
        x = time,
        y = incidence / max_inc
      ),
      linewidth = 1,
      colour = "grey40"
    ) +
    ggplot2::scale_colour_manual(values = state_cols, name = NULL) +
    ggplot2::scale_y_continuous(
      name = "proportion (S / I / R)",
      limits = c(0,1),
      sec.axis = ggplot2::sec_axis(
        ~ . * max_inc,
        name = "incidence (count)"
      )
    ) +
    ggplot2::labs(
      x = "day",
      title = "S, I, R + Daily incidence"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title.y.right = ggplot2::element_text(angle = 90)
    )
  
  # ---------------------------
  # return based on `which`
  # ---------------------------
  
  if (which == "sir") {
    return(p_sir)
  }
  
  if (which == "incidence") {
    return(p_inc)
  }
  
  if (which == "both_side") {
    p_both <- patchwork::wrap_plots(
      p_sir,
      p_inc,
      ncol = 2
    )
    return(p_both)
  }
  
  
  if (which == "overlay") {
    return(p_overlay)
  }
  
  # should never get here
  stop("Internal error: unknown 'which' mode")
}
