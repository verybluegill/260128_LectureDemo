# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

# MSEループ（PDCA）：OM → Obs → Assess → HCR → 次年Catch
mse_loop <- function(TTYear, r, K, B1, q, sigma_CPUE,
                     CR_init = 0.3, management_start = 10,
                     hcr_type = c("hockey", "2k"),
                     hcr_par = list(),
                     under_report_rho = 1.0,
                     fallback_catch = "previous",
                     guard = 1e-8,
                     use_cpp = FALSE) {
  hcr_type <- match.arg(hcr_type)

  # HCRのデフォルト値
  if (hcr_type == "hockey") {
    hcr_par <- modifyList(list(beta = 0.8, Bban_ratio = 0.2, Blim_ratio = 0.5, K_ref_type = "estimated"), hcr_par)
  } else {
    hcr_par <- modifyList(list(n = 5, BT = 0.8, PL = 0.7, PB = 0.0,
                               beta_2k = 0.5, delta = 0.4, lambda = 0.4, clamp = 50), hcr_par)
  }

  # 保存用ベクトル
  B_true <- rep(NA_real_, TTYear)
  C_true <- rep(NA_real_, TTYear)
  C_obs <- rep(NA_real_, TTYear)
  I_true <- rep(NA_real_, TTYear)
  I_obs <- rep(NA_real_, TTYear)
  B_hat <- rep(NA_real_, TTYear)
  K_hat <- rep(NA_real_, TTYear)
  fit_success <- rep(FALSE, TTYear)

  # 初期値
  B_true[1] <- B1
  C_true[1] <- CR_init * B_true[1]

  for (t in 1:(TTYear - 1)) {
    # -----------------------------
    # 1) OM（真の資源動態）
    # -----------------------------
    B_true[t + 1] <- om_step_schaefer(B_true[t], C_true[t], r, K, guard = guard)

    # -----------------------------
    # 2) Obs（CPUE観測）
    # -----------------------------
    obs <- obs_cpue_lognormal(B_true[t], q, sigma_CPUE)
    I_true[t] <- obs$I_true
    I_obs[t] <- obs$I_obs

    # 漁獲の過小報告（rho）
    C_obs[t] <- under_report_rho * C_true[t]

    # -----------------------------
    # 3) Assess + 4) HCR
    # -----------------------------
    if (t >= management_start) {
      # 推定：C++版が使えるなら高速化（Rcppが無ければ自動でR版にフォールバック）
      if (use_cpp) {
        assess <- fit_schaefer_ml_cpp(C_obs[1:t], I_obs[1:t])
      } else {
        assess <- fit_schaefer_ml(C_obs[1:t], I_obs[1:t])
      }

      if (isTRUE(assess$success)) {
        fit_success[t] <- TRUE
        B_hat[t] <- assess$B_hat[t]
        K_hat[t] <- assess$par["K"]
      }

      if (hcr_type == "hockey") {
        K_ref <- if (hcr_par$K_ref_type == "estimated" && is.finite(K_hat[t])) K_hat[t] else K
        Bban <- hcr_par$Bban_ratio * K_ref
        Blim <- hcr_par$Blim_ratio * K_ref

        h_out <- hcr_hockey_stick(B_hat[t], beta = hcr_par$beta, Bban = Bban, Blim = Blim)
        C_next <- h_out$C_next
      } else {
        h_out <- hcr_new2kei_2k(C_obs[1:t], I_obs[1:t],
                                n = hcr_par$n, BT = hcr_par$BT, PL = hcr_par$PL, PB = hcr_par$PB,
                                beta_2k = hcr_par$beta_2k, delta = hcr_par$delta, lambda = hcr_par$lambda,
                                clamp = hcr_par$clamp)
        C_next <- h_out$ABC_t
      }

      # 推定失敗などのフォールバック
      if (!is.finite(C_next) || C_next < 0) {
        if (fallback_catch == "cr") {
          C_next <- CR_init * B_true[t + 1]
        } else {
          C_next <- C_true[t]
        }
      }

      C_true[t + 1] <- C_next
    } else {
      # 管理開始前は決定論C=CR*B
      C_true[t + 1] <- CR_init * B_true[t + 1]
    }
  }

  # 最終年の観測を作成
  obs_last <- obs_cpue_lognormal(B_true[TTYear], q, sigma_CPUE)
  I_true[TTYear] <- obs_last$I_true
  I_obs[TTYear] <- obs_last$I_obs
  C_obs[TTYear] <- under_report_rho * C_true[TTYear]

  list(
    B_true = B_true,
    C_true = C_true,
    C_obs = C_obs,
    I_true = I_true,
    I_obs = I_obs,
    B_hat = B_hat,
    K_hat = K_hat,
    fit_success = fit_success
  )
}
