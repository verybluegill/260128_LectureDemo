# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

# Schaefer型の決定論OM（shape=1）
om_schaefer_sim <- function(TTYear, r, K, B1, CR = 0.3, C_vec = NULL, CR_mult = NULL, guard = 1e-8) {
  # 引数の簡易チェック
  if (TTYear < 2) stop("TTYear must be >= 2")

  if (is.null(C_vec)) {
    C_vec <- rep(NA_real_, TTYear)
  } else {
    if (length(C_vec) != TTYear) stop("length(C_vec) must match TTYear")
  }

  # 年ごとの漁獲係数の揺らぎ（CR_mult=1で決定論）
  if (is.null(CR_mult)) {
    CR_mult <- rep(1, TTYear)
  } else {
    if (length(CR_mult) != TTYear) stop("length(CR_mult) must match TTYear")
  }

  B_true <- rep(NA_real_, TTYear)
  C_true <- rep(NA_real_, TTYear)

  # 初期資源量
  B_true[1] <- B1

  # 年次更新
  for (t in 1:(TTYear - 1)) {
    # ---- 漁獲の決定（決定論） ----
    if (is.na(C_vec[t])) {
      C_true[t] <- CR * CR_mult[t] * B_true[t]
    } else {
      C_true[t] <- C_vec[t]
    }

    # ---- 資源動態（Schaefer：shape=1） ----
    # ここが B_{t+1} = B_t + r*B_t*(1 - B_t/K) - C_t
    B_next <- B_true[t] + r * B_true[t] * (1 - B_true[t] / K) - C_true[t]

    # 数値ガード（負にならないように下限を置く）
    B_true[t + 1] <- max(B_next, guard)
  }

  # 最終年の漁獲も記録（B_Tから決定論で計算）
  if (is.na(C_vec[TTYear])) {
    C_true[TTYear] <- CR * CR_mult[TTYear] * B_true[TTYear]
  } else {
    C_true[TTYear] <- C_vec[TTYear]
  }

  list(B_true = B_true, C_true = C_true)
}

# 1ステップ更新（MSE用）
om_step_schaefer <- function(B_t, C_t, r, K, guard = 1e-8) {
  # ここが B_{t+1} = B_t + r*B_t*(1 - B_t/K) - C_t
  B_next <- B_t + r * B_t * (1 - B_t / K) - C_t
  max(B_next, guard)
}

# Schaeferモデルの予測B系列（推定用）
predict_B_schaefer <- function(C_obs, r, K, B1, guard = 1e-8) {
  TTYear <- length(C_obs)
  B <- rep(NA_real_, TTYear)
  B[1] <- B1

  if (TTYear >= 2) {
    for (t in 1:(TTYear - 1)) {
      B[t + 1] <- om_step_schaefer(B[t], C_obs[t], r, K, guard = guard)
    }
  }

  B
}
