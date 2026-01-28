# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

# ---- 指標計算 ----
calc_rmse <- function(x, y) {
  sqrt(mean((x - y)^2, na.rm = TRUE))
}

calc_relerr_final <- function(B_hat, B_true) {
  (B_hat[length(B_hat)] - B_true[length(B_true)]) / B_true[length(B_true)]
}

calc_aav <- function(x) {
  if (length(x) < 2) return(NA_real_)
  denom <- x[-1] + x[-length(x)]
  numer <- 2 * abs(diff(x))
  frac <- ifelse(denom > 0, numer / denom, 0)
  mean(frac, na.rm = TRUE)
}

calc_mse_metrics <- function(B_true, C_true, K_true, management_start = 10, collapse_ratio = 0.2) {
  idx <- management_start:length(B_true)
  data.frame(
    MeanCatch = mean(C_true[idx], na.rm = TRUE),
    CollapseProb = mean(B_true[idx] < collapse_ratio * K_true, na.rm = TRUE),
    CatchAAV = calc_aav(C_true[idx])
  )
}

# ---- 保存・描画ユーティリティ ----
save_png <- function(path, width = 1000, height = 800, res = 120, plot_fun) {
  png(path, width = width, height = height, res = res)
  on.exit(dev.off())
  plot_fun()
}

# y軸の最小値を0に固定
scale_y_zero <- function() {
  ggplot2::scale_y_continuous(limits = c(0, NA), expand = c(0, 0.05))
}

# ggplotを縦に並べて描画
draw_gg_stack <- function(plots, heights = NULL) {
  n <- length(plots)
  if (n == 0) return(invisible(NULL))
  if (is.null(heights)) heights <- rep(1, n)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, 1, heights = heights)))
  for (i in seq_len(n)) {
    print(plots[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  invisible(NULL)
}

plot_theme_basic <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_ts_case <- function(years, B_true, C_true, I_obs, vline_year = NULL, main_title = "") {
  df_bio <- data.frame(
    Year = rep(years, 2),
    value = c(B_true, I_obs * 20),
    series = rep(c("Biomass", "20 x Index"), each = length(years))
  )
  p_bio <- ggplot2::ggplot(df_bio, ggplot2::aes(x = Year, y = value, color = series)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.6) +
    scale_y_zero() +
    ggplot2::labs(x = "Year", y = "Biomass", title = main_title, color = NULL) +
    ggplot2::scale_color_manual(values = c("Biomass" = "black", "20 x Index" = "darkorange")) +
    plot_theme_basic()

  if (!is.null(vline_year)) {
    p_bio <- p_bio + ggplot2::geom_vline(xintercept = vline_year, linetype = 2, color = "gray40")
  }

  df_catch <- data.frame(Year = years, Catch = C_true)
  p_catch <- ggplot2::ggplot(df_catch, ggplot2::aes(x = Year, y = Catch)) +
    ggplot2::geom_line(color = "black", linewidth = 1) +
    ggplot2::geom_point(color = "black", size = 1.6) +
    scale_y_zero() +
    ggplot2::labs(x = "Year", y = "Catch") +
    plot_theme_basic()

  if (!is.null(vline_year)) {
    p_catch <- p_catch + ggplot2::geom_vline(xintercept = vline_year, linetype = 2, color = "gray40")
  }

  df_index <- data.frame(Year = years, Index = I_obs)
  p_index <- ggplot2::ggplot(df_index, ggplot2::aes(x = Year, y = Index)) +
    ggplot2::geom_point(color = "black", size = 1.6) +
    scale_y_zero() +
    ggplot2::labs(x = "Year", y = "Index") +
    plot_theme_basic()

  if (!is.null(vline_year)) {
    p_index <- p_index + ggplot2::geom_vline(xintercept = vline_year, linetype = 2, color = "gray40")
  }

  draw_gg_stack(list(p_bio, p_catch, p_index))
}

plot_true_vs_hat <- function(years, B_true, B_hat, main_title = "") {
  df <- data.frame(
    Year = rep(years, 2),
    Biomass = c(B_true, B_hat),
    series = rep(c("True", "Estimated"), each = length(years))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Year, y = Biomass, color = series)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.6) +
    scale_y_zero() +
    ggplot2::labs(x = "Year", y = "Biomass", title = main_title, color = NULL) +
    ggplot2::scale_color_manual(values = c("True" = "black", "Estimated" = "red")) +
    plot_theme_basic()

  print(p)
  invisible(p)
}

plot_sensitivity_heatmap <- function(summary_df) {
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = TTYear, y = sigma_CPUE, fill = RMSE_B_mean)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(RMSE_B_mean, 2)), size = 3) +
    ggplot2::scale_x_continuous(breaks = sort(unique(summary_df$TTYear))) +
    ggplot2::scale_y_continuous(breaks = sort(unique(summary_df$sigma_CPUE)), limits = c(0, NA), expand = c(0, 0.05)) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::labs(x = "TTYear", y = "sigma_Index", fill = "RMSE") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}

plot_hockey_shape <- function(K_ref, beta = 0.8, Bban_ratio = 0.2, Blim_ratio = 0.5) {
  Bban <- Bban_ratio * K_ref
  Blim <- Blim_ratio * K_ref
  B_seq <- seq(0, K_ref * 1.2, length.out = 200)

  df <- hcr_hockey_shape(B_seq, beta = beta, Bban = Bban, Blim = Blim)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = B, y = h)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = c(Bban, Blim), linetype = 2, color = "gray40") +
    scale_y_zero() +
    ggplot2::labs(x = "Biomass", y = "Harvest rate", title = "Hockey-stick HCR") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}

plot_2k_shape <- function(AAV_t = 0.2, BT = 0.8, PL = 0.7, PB = 0.0,
                          beta_2k = 0.5, delta = 0.4, lambda = 0.4) {
  D_seq <- seq(0, 1, length.out = 200)
  ABC_ratio <- sapply(D_seq, function(d) {
    hcr_new2kei_2k_from_D(d, AAV_t = AAV_t, Cbar_t = 1,
                          BT = BT, PL = PL, PB = PB,
                          beta_2k = beta_2k, delta = delta, lambda = lambda)$ABC_ratio
  })

  df <- data.frame(D = D_seq, ABC_ratio = ABC_ratio)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = D, y = ABC_ratio)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = c(PL * BT, PB * BT), linetype = 2, color = "gray40") +
    scale_y_zero() +
    ggplot2::labs(x = "D (standardized Index)", y = "ABC / Cbar", title = "Type2 rule") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}
