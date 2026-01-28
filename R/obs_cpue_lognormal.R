# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

# CPUE観測：対数正規（バイアス補正あり）
obs_cpue_lognormal <- function(B_true, q, sigma_CPUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  I_true <- q * B_true
  # e_t ~ N(-0.5*sigma^2, sigma)
  e_t <- rnorm(length(B_true), mean = -0.5 * sigma_CPUE^2, sd = sigma_CPUE)
  I_obs <- I_true * exp(e_t)

  list(I_true = I_true, I_obs = I_obs, e_t = e_t)
}

# 検証用：固定Bで平均がq*Bに一致するか確認
check_cpue_lognormal_mean <- function(B_fixed, q, sigma_CPUE, nsim = 100000, seed = 123) {
  set.seed(seed)

  I_true <- q * B_fixed
  e_t <- rnorm(nsim, mean = -0.5 * sigma_CPUE^2, sd = sigma_CPUE)
  I_obs <- I_true * exp(e_t)
  ratio <- mean(I_obs) / I_true

  list(ratio = ratio, I_obs = I_obs)
}