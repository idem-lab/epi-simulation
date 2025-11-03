plot_sirs <- function(
    sim,
    which = c("both_side","overlay","sir","incidence"),
    per_million = FALSE
) {
  which <- match.arg(which)
  
  # ---- basic checks for required series ----
  need <- c("time","S","I","R","incidence")
  miss <- setdiff(need, names(sim))
  if (length(miss)) stop("sim is missing: ", paste(miss, collapse=", "))
  
  # ---- get pop from sim$params$pop if needed ----
  if (per_million) {
    pop_ok <-
      !is.null(sim$params$pop) &&
      is.numeric(sim$params$pop) &&
      length(sim$params$pop) == 1 &&
      is.finite(sim$params$pop) &&
      sim$params$pop > 0
    
    if (!pop_ok) {
      stop(
        "per_million = TRUE but sim$params$pop is not a single positive number.\n",
        "Please make sure sim$params$pop exists and is > 0."
      )
    }
    
    pop_val <- sim$params$pop
  } else {
    pop_val <- NA_real_  # won't be used
  }
  
  # ---------------------------------
  # prepare data, scale incidence if requested
  # ---------------------------------
  if (per_million) {
    incidence_val <- (sim$incidence / pop_val) * 1e6
    incidence_label <- "new infections"
  } else {
    incidence_val <- sim$incidence
    incidence_label <- "new infections (count)"
  }
  
  df_sir <- data.frame(
    time = sim$time,
    S    = sim$S,
    I    = sim$I,
    R    = sim$R
  )
  
  df_inc <- data.frame(
    time      = sim$time,
    incidence = incidence_val
  )
  
  # long format for S/I/R
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
  # p_sir (S, I, R proportions)
  # ---------------------------
  p_sir <- ggplot2::ggplot(
    df_long,
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
  
  # ---------------------------
  # p_inc (incidence panel)
  # ---------------------------
  p_inc <- ggplot2::ggplot(
    df_inc,
    ggplot2::aes(x = time, y = incidence)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "day",
      y = incidence_label,
      title = "Daily incidence"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  # ---------------------------
  # p_overlay (dual axis)
  # ---------------------------
  # For overlay: rescale incidence to [0,1] so it sits on same panel
  # as S/I/R proportions. Then expose a secondary axis that maps back.
  max_inc <- max(df_inc$incidence, na.rm = TRUE)
  if (!is.finite(max_inc) || max_inc == 0) {
    max_inc <- 1
  }
  
  p_overlay <- ggplot2::ggplot() +
    # S / I / R
    ggplot2::geom_line(
      data = df_long,
      ggplot2::aes(x = time, y = value, colour = state),
      linewidth = 1
    ) +
    # scaled incidence
    ggplot2::geom_line(
      data = df_inc,
      ggplot2::aes(
        x = time,
        y = incidence / max_inc
      ),
      linewidth = 1,
      colour = "grey60"
    ) +
    ggplot2::scale_colour_manual(values = state_cols, name = NULL) +
    ggplot2::scale_y_continuous(
      name = "proportion (S / I / R)",
      limits = c(0,1),
      sec.axis = ggplot2::sec_axis(
        ~ . * max_inc,
        name = incidence_label
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
  # return view
  # ---------------------------
  if (which == "sir") {
    return(p_sir)
  }
  
  if (which == "incidence") {
    return(p_inc)
  }
  
  if (which == "both_side") {
    # use patchwork to arrange without relying on +
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
  
  stop("Internal error: unknown 'which' mode")
}
