# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

# 新2系（2K）ルール
hcr_new2kei_2k <- function(C_obs, CPUE_obs, n = 5, BT = 0.8, PL = 0.7, PB = 0.0,
                           beta_2k = 0.5, delta = 0.4, lambda = 0.4, clamp = 50) {
  t <- length(CPUE_obs)
  if (t == 0) stop("CPUE_obs length is 0")

  # 過去n年平均漁獲
  idx_n <- max(1, t - n + 1):t
  Cbar_t <- mean(C_obs[idx_n], na.rm = TRUE)

  # CPUE標準化
  mu_t <- mean(CPUE_obs, na.rm = TRUE)
  sd_t <- sd(CPUE_obs, na.rm = TRUE)
  if (!is.finite(sd_t) || sd_t <= 0) {
    D_t <- 0.5
  } else {
    z_t <- (CPUE_obs[t] - mu_t) / sd_t
    D_t <- pnorm(z_t)
  }

  # AAV（CPUEの年変動指標）
  if (t < 2) {
    AAV_t <- 0
  } else {
    denom <- CPUE_obs[2:t] + CPUE_obs[1:(t - 1)]
    numer <- 2 * abs(CPUE_obs[2:t] - CPUE_obs[1:(t - 1)])
    frac <- ifelse(denom > 0, numer / denom, 0)
    AAV_t <- mean(frac, na.rm = TRUE)
  }

  # 閾値
  BL <- PL * BT
  BB <- PB * BT

  # k_t の決定
  if (D_t > BL) {
    k_t <- beta_2k
  } else if (D_t > BB && D_t <= BL) {
    k_t <- beta_2k + delta * exp(lambda * log(AAV_t^2 + 1)) * (BL - D_t) / (D_t - BB)
  } else {
    k_t <- Inf
  }

  # ABCの計算
  if (!is.finite(Cbar_t) || Cbar_t < 0) {
    ABC_t <- NA_real_
  } else if (is.infinite(k_t)) {
    ABC_t <- 0
  } else {
    expo <- k_t * (D_t - BT)
    expo <- max(min(expo, clamp), -clamp)
    ABC_t <- exp(expo) * Cbar_t
  }

  ABC_ratio <- ifelse(Cbar_t > 0, ABC_t / Cbar_t, NA_real_)

  list(
    D_t = D_t,
    AAV_t = AAV_t,
    Cbar_t = Cbar_t,
    k_t = k_t,
    ABC_t = ABC_t,
    ABC_ratio = ABC_ratio,
    BL = BL,
    BB = BB
  )
}

# 形状図用：D_t を直接入れて ABC/Cbar を計算する簡易版
hcr_new2kei_2k_from_D <- function(D_t, AAV_t = 0.2, Cbar_t = 1,
                                  BT = 0.8, PL = 0.7, PB = 0.0,
                                  beta_2k = 0.5, delta = 0.4, lambda = 0.4, clamp = 50) {
  BL <- PL * BT
  BB <- PB * BT

  if (D_t > BL) {
    k_t <- beta_2k
  } else if (D_t > BB && D_t <= BL) {
    k_t <- beta_2k + delta * exp(lambda * log(AAV_t^2 + 1)) * (BL - D_t) / (D_t - BB)
  } else {
    k_t <- Inf
  }

  if (is.infinite(k_t)) {
    ABC_t <- 0
  } else {
    expo <- k_t * (D_t - BT)
    expo <- max(min(expo, clamp), -clamp)
    ABC_t <- exp(expo) * Cbar_t
  }

  ABC_ratio <- ifelse(Cbar_t > 0, ABC_t / Cbar_t, NA_real_)

  list(ABC_ratio = ABC_ratio, k_t = k_t, BL = BL, BB = BB)
}