
mse_loop <- function(TTYear, r, K, B1, q, sigma_CPUE,
                     CR_init = 0.3, management_start = 10,
                     hcr_type = c("hockey", "2k"),
                     hcr_par = list(),
                     B_input = NULL,
                     under_report_rho = 1.0,
                     fallback_catch = "previous",
                     guard = 1e-8,
                     biomass_floor_ratio = 0.05,
                     use_cpp = FALSE,
                     brake = list(u_max = Inf, g_max = Inf),
                     debug = FALSE) {
  hcr_type <- match.arg(hcr_type)

  if (hcr_type == "hockey") {
    hcr_par <- modifyList(list(beta = 0.8, Bban_ratio = 0.2, Blim_ratio = 0.5, K_ref_type = "estimated"), hcr_par)
  } else {
    hcr_par <- modifyList(list(n = 5, BT = 0.8, PL = 0.7, PB = 0.0,
                               beta_2k = 0.5, delta = 0.4, lambda = 0.4, clamp = 50), hcr_par)
  }

  if (is.null(brake$u_max)) brake$u_max <- Inf
  if (is.null(brake$g_max)) brake$g_max <- Inf

  B_true <- rep(NA_real_, TTYear)
  C_true <- rep(NA_real_, TTYear)
  C_obs <- rep(NA_real_, TTYear)
  I_true <- rep(NA_real_, TTYear)
  I_obs <- rep(NA_real_, TTYear)
  B_hat <- rep(NA_real_, TTYear)
  K_hat <- rep(NA_real_, TTYear)
  fit_success <- rep(FALSE, TTYear)

  B_true[1] <- B1
  C_true[1] <- CR_init * B_true[1]
  B_floor <- biomass_floor_ratio * K

  for (t in 1:(TTYear - 1)) {
    B_true[t + 1] <- om_step_schaefer(B_true[t], C_true[t], r, K, guard = guard)

    obs <- obs_cpue_lognormal(B_true[t], q, sigma_CPUE)
    I_true[t] <- obs$I_true
    I_obs[t] <- obs$I_obs

    C_obs[t] <- under_report_rho * C_true[t]

    if (t >= management_start) {
      if (!is.null(B_input) && length(B_input) >= t && is.finite(B_input[t])) {
        assess <- list(success = FALSE, B_hat = rep(NA_real_, t), par = c(K = NA_real_, r = NA_real_))
        fit_success[t] <- FALSE
        B_hat[t] <- NA_real_
        K_hat[t] <- NA_real_
        B_for_hcr <- B_input[t]
        B_source <- "input"
      } else {
      if (use_cpp) {
        assess <- tryCatch(
          fit_schaefer_ml_cpp(C_obs[1:t], I_obs[1:t]),
          error = function(e) fit_schaefer_ml(C_obs[1:t], I_obs[1:t])
        )
      } else {
        assess <- fit_schaefer_ml(C_obs[1:t], I_obs[1:t])
      }

      b_hat_t <- if (!is.null(assess$B_hat) && length(assess$B_hat) >= t) assess$B_hat[t] else NA_real_
      assess_ok <- isTRUE(assess$success) && is.finite(b_hat_t)
      if (assess_ok) {
        fit_success[t] <- TRUE
        B_hat[t] <- b_hat_t
        K_hat[t] <- assess$par["K"]
      } else {
        fit_success[t] <- FALSE
        B_hat[t] <- NA_real_
      }

      B_for_hcr <- if (is.finite(B_hat[t])) B_hat[t] else B_true[t]
      B_source <- if (is.finite(B_hat[t])) "hat" else "true"
      }

      if (hcr_type == "hockey") {
        K_ref <- if (hcr_par$K_ref_type == "estimated" && is.finite(K_hat[t])) K_hat[t] else K
        Bban <- hcr_par$Bban_ratio * K_ref
        Blim <- hcr_par$Blim_ratio * K_ref
        r_hat_t <- if (!is.null(assess$par) && "r" %in% names(assess$par) && is.finite(assess$par["r"])) assess$par["r"] else NA_real_
        Fmsy_ref <- if (is.finite(r_hat_t)) r_hat_t / 2 else r / 2

        h_out <- hcr_hockey_stick(B_for_hcr, beta = hcr_par$beta, Bban = Bban, Blim = Blim, Fmsy = Fmsy_ref)
        C_next <- h_out$C_next
      } else {
        h_out <- hcr_new2kei_2k(C_obs[1:t], I_obs[1:t],
                                n = hcr_par$n, BT = hcr_par$BT, PL = hcr_par$PL, PB = hcr_par$PB,
                                beta_2k = hcr_par$beta_2k, delta = hcr_par$delta, lambda = hcr_par$lambda,
                                clamp = hcr_par$clamp)
        C_next <- h_out$ABC_t
      }

      C_next_raw <- C_next

      if (!is.finite(C_next) || C_next < 0) {
        if (fallback_catch == "cr") {
          C_next <- CR_init * B_true[t + 1]
        } else {
          C_next <- C_true[t]
        }
      }

      C_next_pre_brake <- C_next
      clipped_u <- FALSE
      clipped_g <- FALSE
      if (is.finite(C_next) && is.finite(B_for_hcr) && is.finite(brake$u_max)) {
        C_u_max <- brake$u_max * B_for_hcr
        if (C_next > C_u_max) {
          C_next <- C_u_max
          clipped_u <- TRUE
        }
      }
      if (is.finite(C_next) && is.finite(C_true[t]) && is.finite(brake$g_max)) {
        C_g_max <- (1 + brake$g_max) * C_true[t]
        if (C_next > C_g_max) {
          C_next <- C_g_max
          clipped_g <- TRUE
        }
      }
      C_next_post_brake <- C_next
      clipped_brake <- (isTRUE(clipped_u) || isTRUE(clipped_g))

      C_max <- max(0, B_true[t] + r * B_true[t] * (1 - B_true[t] / K) - guard)
      C_next_pre_clip <- C_next
      C_next <- min(max(C_next, 0), C_max)

      prod_t <- r * B_true[t] * (1 - B_true[t] / K)
      C_floor_max <- max(0, B_true[t] + prod_t - B_floor)
      C_next_pre_floor <- C_next
      C_next <- min(C_next, C_floor_max)
      clipped_floor <- (C_next_pre_floor != C_next)

      if (isTRUE(debug)) {
        cat(sprintf(
          "t=%d, success=%s, B_source=%s, B_hat=%s, B_true=%s, C_next_raw=%s, C_next_pre_brake=%s, C_next_post_brake=%s, clipped_u=%s, clipped_g=%s, clipped_brake=%s, C_next_pre_clip=%s, C_next=%s, C_max=%s, B_floor=%s, C_floor_max=%s, clipped_floor=%s\n",
          t,
          if (isTRUE(assess$success)) "TRUE" else "FALSE",
          B_source,
          ifelse(is.finite(B_hat[t]), format(B_hat[t], digits = 6), as.character(B_hat[t])),
          format(B_true[t], digits = 6),
          ifelse(is.finite(C_next_raw), format(C_next_raw, digits = 6), as.character(C_next_raw)),
          ifelse(is.finite(C_next_pre_brake), format(C_next_pre_brake, digits = 6), as.character(C_next_pre_brake)),
          ifelse(is.finite(C_next_post_brake), format(C_next_post_brake, digits = 6), as.character(C_next_post_brake)),
          if (isTRUE(clipped_u)) "TRUE" else "FALSE",
          if (isTRUE(clipped_g)) "TRUE" else "FALSE",
          if (isTRUE(clipped_brake)) "TRUE" else "FALSE",
          ifelse(is.finite(C_next_pre_clip), format(C_next_pre_clip, digits = 6), as.character(C_next_pre_clip)),
          format(C_next, digits = 6),
          format(C_max, digits = 6),
          format(B_floor, digits = 6),
          format(C_floor_max, digits = 6),
          if (isTRUE(clipped_floor)) "TRUE" else "FALSE"
        ))
      }

      C_true[t + 1] <- C_next
    } else {
      C_true[t + 1] <- CR_init * B_true[t + 1]
    }
  }

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
