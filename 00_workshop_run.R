# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

############################################################
# 00_workshop_run.R
# 蛻晏ｿ・・髄縺托ｼ夊ｳ・ｺ占ｩ穂ｾ｡繝ｻHCR繝ｻMSE貍皮ｿ抵ｼ井ｸ翫°繧蛾・↓螳溯｡鯉ｼ・############################################################

# ==========================================================
# 1) 貅門ｙ・壹ヱ繝・こ繝ｼ繧ｸ縲《eed縲《ource縲∝・蜉帙ヵ繧ｩ繝ｫ繝
# ==========================================================
rm(list = ls())
set.seed(123)
options(stringsAsFactors = FALSE)

# 蜿ら・繧ｳ繝ｼ繝峨・隱ｭ縺ｿ霎ｼ縺ｿ・亥､画焚蜷阪・髮ｰ蝗ｲ豌励ｒ邯ｭ謖√☆繧九◆繧・ｼ・ref_code_path <- "C:/Users/rmits/Documents/1_Desk/260128_LectureDemo/CreatingTD.Rmd"
ref_lines <- readLines(ref_code_path, warn = FALSE)
cat("蜿ら・繧ｳ繝ｼ繝芽｡梧焚:", length(ref_lines), "\n")
# 縺薙％縺ｧ菴ｿ繧上ｌ縺ｦ縺・ｋ螟画焚蜷搾ｼ・, Catch, cpue_obs 縺ｪ縺ｩ・峨ｒ譛ｬ謨呎攝縺ｧ繧りｸ剰･ｲ

# 髢｢謨ｰ縺ｮ隱ｭ縺ｿ霎ｼ縺ｿ
source("R/om_schaefer.R")
source("R/obs_cpue_lognormal.R")
source("R/fit_schaefer_ml.R")
source("R/hcr_hockey_stick.R")
source("R/hcr_new2kei_2k.R")
source("R/mse_loop.R")
source("R/metrics_plots.R")

# 蜃ｺ蜉帙ヵ繧ｩ繝ｫ繝
out_fig <- file.path("output", "fig")
out_tab <- file.path("output", "tab")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

# ==========================================================
# 2) CPUE縺悟ｯｾ謨ｰ豁｣隕上↓縺ｪ縺｣縺ｦ縺・ｋ縺薙→縺ｮ讀懆ｨｼ
# ==========================================================
# 逶ｮ逧・ｼ哘[I_obs | B] = q*B 繧堤｢ｺ隱・B_fixed <- 500
q_check <- 0.05
sigma_check <- 0.2

chk <- check_cpue_lognormal_mean(B_fixed, q_check, sigma_check, nsim = 50000, seed = 999)
cat("CPUE check: mean(I_obs) / (q*B) = ", round(chk$ratio, 4), "\n")

# 蝗ｳ・唔_obs / (q*B) 縺ｮ蛻・ｸ・save_png(file.path(out_fig, "cpue_lognormal_check.png"), width = 900, height = 600, plot_fun = function() {
  ratio_vec <- chk$I_obs / (q_check * B_fixed)
  hist(ratio_vec, breaks = 50, main = "Lognormal CPUE check",
       xlab = "I_obs / (q*B)", col = "gray80", border = "gray40")
  abline(v = 1, col = "red", lwd = 2)
})

# ==========================================================
# 3) 1繧ｱ繝ｼ繧ｹ縺ｮ繝・・繧ｿ逕滓・・井ｾ具ｼ啜TYear=50, sigma=0.1・・# ==========================================================
# 繝代Λ繝｡繝ｼ繧ｿ險ｭ螳・TTYear <- 50
TTYear <- 50
r <- 0.6
K <- 1000
D1 <- 0.7
B1 <- D1 * K
CR <- 0.2
q <- 0.05
sigma_CPUE <- 0.1

# 高速化の設定（Rcppが無い場合は自動でR版にフォールバック）
use_cpp <- TRUE
if (use_cpp && !requireNamespace("Rcpp", quietly = TRUE)) {
  message("Rcpp が見つからないため、use_cpp = FALSE にします。")
  use_cpp <- FALSE
}

# OM・・chaefer・峨〒逵溘・雉・ｺ宣㍼縺ｨ貍∫佐繧堤函謌・om <- om_schaefer_sim(TTYear = TTYear, r = r, K = K, B1 = B1, CR = CR)
B_true <- om$B_true
Catch <- om$C_true

# 蜿り・さ繝ｼ繝峨〒縺ｯ貍∫佐縺ｫ荵ｱ謨ｰ繧貞・繧後※縺・◆縺後∵蕗譚舌〒縺ｯ豎ｺ螳夊ｫ悶↓蝗ｺ螳・# Catch[t] <- CR * B[t] * runif(1, 0.7, 1.5)  # 竊・蜿り・さ繝ｼ繝峨・謠ｺ繧峨℃
# Catch[t] <- CR * B[t]                       # 竊・譛ｬ謨呎攝縺ｯ縺薙ｌ繧呈治逕ｨ

# CPUE隕ｳ貂ｬ・亥ｯｾ謨ｰ豁｣隕上√ヰ繧､繧｢繧ｹ陬懈ｭ｣縺ゅｊ・・obs <- obs_cpue_lognormal(B_true, q, sigma_CPUE)
cpue_obs <- obs$I_obs

Year <- 1975:(1975 + TTYear - 1)

# 繝・・繧ｿ縺ｾ縺ｨ繧・case_data <- data.frame(
  Year = Year,
  Biomass = B_true,
  Catch = Catch,
  CPUE = cpue_obs
)

# 蝗ｳ・咤, Catch, CPUE
save_png(file.path(out_fig, "case_timeseries.png"), width = 900, height = 900, plot_fun = function() {
  plot_ts_case(Year, B_true, Catch, cpue_obs, main_title = "Single-case OM + Obs")
})

# ==========================================================
# 4) Schaefer謗ｨ螳夲ｼ域怙蟆､・・# ==========================================================
# 隕ｳ貂ｬ貍∫佐・磯℃蟆丞ｱ蜻翫↑縺励・蝓ｺ譛ｬ繧ｱ繝ｼ繧ｹ・・Catch_obs <- Catch

fit <- if (use_cpp) fit_schaefer_ml_cpp(Catch_obs, cpue_obs) else fit_schaefer_ml(Catch_obs, cpue_obs)

# 謗ｨ螳壹・謌仙凄
cat("Fit success:", fit$success, "  converged:", fit$converged, "\n")

# 逵溷､縺ｨ謗ｨ螳壹・豈碑ｼ・ｼ・邉ｻ蛻暦ｼ・save_png(file.path(out_fig, "fit_biomass.png"), width = 900, height = 600, plot_fun = function() {
  plot_true_vs_hat(Year, B_true, fit$B_hat, main_title = "True vs Estimated Biomass")
})

# 繝代Λ繝｡繝ｼ繧ｿ陦ｨ・育悄蛟､縺ｨ謗ｨ螳壼､・・fit_table <- data.frame(
  parameter = c("r", "K", "q", "B1", "sigma_I"),
  true = c(r, K, q, B1, sigma_CPUE),
  estimate = c(fit$par["r"], fit$par["K"], fit$par["q"], fit$par["B1"], fit$par["sigma_I"])
)
write.csv(fit_table, file.path(out_tab, "fit_params.csv"), row.names = FALSE)

# ==========================================================
# 5) 諢溷ｺｦ隗｣譫撰ｼ・TYear ﾃ・sigma_CPUE・・# ==========================================================
TTYear_vec <- c(30, 50, 100)
sigma_vec <- c(0.05, 0.10, 0.20, 0.30, 0.50)
nrep <- 20  # 險育ｮ苓ｲ闕ｷ繧剃ｸ九￡繧九◆繧∝ｰ代↑繧・ｼ亥ｿ・ｦ√↑繧・50/100 縺ｫ謌ｻ縺呻ｼ・
sens_list <- list()
idx <- 1

for (TT in TTYear_vec) {
  # 豎ｺ螳夊ｫ悶↑縺ｮ縺ｧB_true縺ｯ蝗ｺ螳夲ｼ・PUE縺縺台ｹｱ謨ｰ・・  om_tmp <- om_schaefer_sim(TTYear = TT, r = r, K = K, B1 = D1 * K, CR = CR)

  for (sg in sigma_vec) {
    for (rep in 1:nrep) {
      obs_tmp <- obs_cpue_lognormal(om_tmp$B_true, q, sg)
      fit_tmp <- if (use_cpp) fit_schaefer_ml_cpp(om_tmp$C_true, obs_tmp$I_obs) else fit_schaefer_ml(om_tmp$C_true, obs_tmp$I_obs)

      if (isTRUE(fit_tmp$success)) {
        rmse_b <- calc_rmse(fit_tmp$B_hat, om_tmp$B_true)
        relerr <- calc_relerr_final(fit_tmp$B_hat, om_tmp$B_true)
        success <- 1
      } else {
        rmse_b <- NA_real_
        relerr <- NA_real_
        success <- 0
      }

      sens_list[[idx]] <- data.frame(
        TTYear = TT,
        sigma_CPUE = sg,
        rep = rep,
        RMSE_B = rmse_b,
        RelErr_final = relerr,
        success = success
      )
      idx <- idx + 1
    }
  }
}

sens_df <- do.call(rbind, sens_list)

# 縺ｾ縺ｨ繧∬｡ｨ
sens_summary <- aggregate(RMSE_B ~ TTYear + sigma_CPUE, data = sens_df, FUN = function(x) mean(x, na.rm = TRUE))
colnames(sens_summary)[3] <- "RMSE_B_mean"

sens_summary_sd <- aggregate(RMSE_B ~ TTYear + sigma_CPUE, data = sens_df, FUN = function(x) sd(x, na.rm = TRUE))
colnames(sens_summary_sd)[3] <- "RMSE_B_sd"

sens_summary_rel <- aggregate(RelErr_final ~ TTYear + sigma_CPUE, data = sens_df, FUN = function(x) mean(x, na.rm = TRUE))
colnames(sens_summary_rel)[3] <- "RelErr_final_mean"

sens_fail <- aggregate(success ~ TTYear + sigma_CPUE, data = sens_df, FUN = function(x) 1 - mean(x, na.rm = TRUE))
colnames(sens_fail)[3] <- "FailRate"

sens_summary <- merge(sens_summary, sens_summary_sd, by = c("TTYear", "sigma_CPUE"))
sens_summary <- merge(sens_summary, sens_summary_rel, by = c("TTYear", "sigma_CPUE"))
sens_summary <- merge(sens_summary, sens_fail, by = c("TTYear", "sigma_CPUE"))

write.csv(sens_summary, file.path(out_tab, "sensitivity_summary.csv"), row.names = FALSE)

# 蝗ｳ・啌MSE縺ｮ繝偵・繝医・繝・・
save_png(file.path(out_fig, "sensitivity_rmse_heatmap.png"), width = 900, height = 600, plot_fun = function() {
  plot_sensitivity_heatmap(sens_summary)
})

# ==========================================================
# 6) HCR縺ｮ蠖｢縺ｮ蝗ｳ遉ｺ・・ockey 縺ｨ 2K・・# ==========================================================
K_ref <- K

save_png(file.path(out_fig, "hcr_hockey_shape.png"), width = 900, height = 600, plot_fun = function() {
  plot_hockey_shape(K_ref = K_ref, beta = 0.8, Bban_ratio = 0.3, Blim_ratio = 0.6)
})

save_png(file.path(out_fig, "hcr_2k_shape.png"), width = 900, height = 600, plot_fun = function() {
  plot_2k_shape(AAV_t = 0.2, BT = 0.75, PL = 0.7, PB = 0.0,
                beta_2k = 0.6, delta = 0.3, lambda = 0.3)
})

# ==========================================================
# 7) MSE螳溯｡鯉ｼ喇ockey vs 2K
# ==========================================================
management_start <- 10

mse_hockey <- mse_loop(
  TTYear = TTYear, r = r, K = K, B1 = B1, q = q, sigma_CPUE = sigma_CPUE,
  CR_init = CR, management_start = management_start,
  hcr_type = "hockey",
  hcr_par = list(beta = 0.8, Bban_ratio = 0.3, Blim_ratio = 0.6, K_ref_type = "fixed"),
  under_report_rho = 1.0,
  use_cpp = use_cpp
)

mse_2k <- mse_loop(
  TTYear = TTYear, r = r, K = K, B1 = B1, q = q, sigma_CPUE = sigma_CPUE,
  CR_init = CR, management_start = management_start,
  hcr_type = "2k",
  hcr_par = list(n = 5, BT = 0.75, PL = 0.7, PB = 0.0,
                 beta_2k = 0.6, delta = 0.3, lambda = 0.3, clamp = 50),
  under_report_rho = 1.0,
  use_cpp = use_cpp
)

# 謖・ｨ呵ｨ育ｮ・mse_sum_h <- calc_mse_metrics(mse_hockey$B_true, mse_hockey$C_true, K_true = K, management_start = management_start)
mse_sum_h$HCR <- "hockey"

mse_sum_2k <- calc_mse_metrics(mse_2k$B_true, mse_2k$C_true, K_true = K, management_start = management_start)
mse_sum_2k$HCR <- "2k"

mse_summary <- rbind(mse_sum_h, mse_sum_2k)
write.csv(mse_summary, file.path(out_tab, "mse_summary.csv"), row.names = FALSE)

# 蝗ｳ・咤, Catch, CPUE 繧呈凾邉ｻ蛻励〒豈碑ｼ・ｼ育ｮ｡逅・幕蟋句ｹｴ縺ｫ邵ｦ邱夲ｼ・save_png(file.path(out_fig, "mse_timeseries_compare.png"), width = 900, height = 900, plot_fun = function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

  plot(Year, mse_hockey$B_true, type = "l", xlab = "Year", ylab = "Biomass", main = "Biomass (MSE)")
  lines(Year, mse_2k$B_true, col = "blue", lwd = 2)
  abline(v = Year[management_start], lty = 2, col = "gray40")
  legend("topright", legend = c("Hockey", "2K"), col = c("black", "blue"), lty = 1, bty = "n")

  plot(Year, mse_hockey$C_true, type = "l", xlab = "Year", ylab = "Catch", main = "Catch (MSE)")
  lines(Year, mse_2k$C_true, col = "blue", lwd = 2)
  abline(v = Year[management_start], lty = 2, col = "gray40")

  plot(Year, mse_hockey$I_obs, type = "p", xlab = "Year", ylab = "CPUE", main = "Observed CPUE (MSE)")
  points(Year, mse_2k$I_obs, col = "blue", pch = 1)
  abline(v = Year[management_start], lty = 2, col = "gray40")
})

# ==========================================================
# 8) ・医が繝励す繝ｧ繝ｳ・韻atch驕主ｰ丞ｱ蜻翫・MSE
# ==========================================================
# TRUE縺ｫ縺吶ｋ縺ｨ縲〉ho<1 縺ｮ豈碑ｼ・ｒ螳溯｡後☆繧・# do_underreport <- TRUE
# rho_vec <- c(0.8, 0.6)

# if (do_underreport) {
#   res_list <- list()
#   k <- 1
#   for (rho in rho_vec) {
#     tmp_h <- mse_loop(TTYear, r, K, B1, q, sigma_CPUE, CR, management_start,
#                       hcr_type = "hockey",
#                       hcr_par = list(beta = 0.8, Bban_ratio = 0.3, Blim_ratio = 0.6, K_ref_type = "fixed"),
#                       under_report_rho = rho,
#                       use_cpp = use_cpp)
#     tmp_2k <- mse_loop(TTYear, r, K, B1, q, sigma_CPUE, CR, management_start,
#                        hcr_type = "2k",
#                        hcr_par = list(n = 5, BT = 0.75, PL = 0.7, PB = 0.0,
#                                       beta_2k = 0.6, delta = 0.3, lambda = 0.3, clamp = 50),
#                        under_report_rho = rho,
#                        use_cpp = use_cpp)
#
#     s1 <- calc_mse_metrics(tmp_h$B_true, tmp_h$C_true, K_true = K, management_start = management_start)
#     s1$HCR <- "hockey"
#     s1$rho <- rho
#
#     s2 <- calc_mse_metrics(tmp_2k$B_true, tmp_2k$C_true, K_true = K, management_start = management_start)
#     s2$HCR <- "2k"
#     s2$rho <- rho
#
#     res_list[[k]] <- s1; k <- k + 1
#     res_list[[k]] <- s2; k <- k + 1
#   }
#
#   under_summary <- do.call(rbind, res_list)
#   write.csv(under_summary, file.path(out_tab, "mse_summary_underreport.csv"), row.names = FALSE)
# }

# ==========================================================
# 9) 逕滓・迚ｩ縺ｮ繝代せ荳隕ｧ
# ==========================================================
cat("\n--- output/fig ---\n")
print(list.files(out_fig, full.names = TRUE))
cat("\n--- output/tab ---\n")
print(list.files(out_tab, full.names = TRUE))

############################################################
# 縺薙％縺ｾ縺ｧ
############################################################
