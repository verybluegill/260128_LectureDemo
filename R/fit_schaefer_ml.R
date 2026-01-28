# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

# Schaeferモデルの最尤推定（決定論、CPUE誤差のみ）
fit_schaefer_ml <- function(C_obs, CPUE_obs, start_list = NULL, control = list(maxit = 1000)) {
  TTYear <- length(C_obs)
  if (length(CPUE_obs) != TTYear) stop("length mismatch: C_obs and CPUE_obs")

  # CPUEの正値チェック（対数正規の前提）
  idx <- which(is.finite(CPUE_obs) & CPUE_obs > 0)
  if (length(idx) < 5) {
    return(list(
      success = FALSE,
      message = "CPUE_obs has too few positive values",
      par = NA,
      B_hat = rep(NA_real_, TTYear),
      nll = NA_real_,
      converged = FALSE
    ))
  }

  # 負の対数尤度
  nll_fun <- function(par) {
    r <- exp(par[1])
    K <- exp(par[2])
    q <- exp(par[3])
    B1 <- exp(par[4])
    sigma_I <- exp(par[5])

    # 予測B系列（決定論）
    B_pred <- predict_B_schaefer(C_obs, r, K, B1)
    if (any(!is.finite(B_pred)) || any(B_pred <= 0)) return(1e12)

    # log(CPUE_obs_t) ~ N(log(q*B_t) - 0.5*sigma^2, sigma^2)
    mu <- log(q * B_pred[idx]) - 0.5 * sigma_I^2
    ll <- dnorm(log(CPUE_obs[idx]), mean = mu, sd = sigma_I, log = TRUE)
    -sum(ll)
  }

  # 初期値（簡易）
  q0 <- 0.05
  first <- idx[1]
  B1_approx <- max(CPUE_obs[first] / q0, 1e-6)  
  K0 <- max(B1_approx * 2, max(C_obs, na.rm = TRUE) * 5, 1)
  r0 <- 0.5
  sigma0 <- max(0.05, sd(log(CPUE_obs[idx]), na.rm = TRUE))
  if (!is.finite(sigma0)) sigma0 <- 0.2

  base_start <- c(log(r0), log(K0), log(q0), log(B1_approx), log(sigma0))

  if (is.null(start_list)) {
    # 単純なマルチスタート（失敗率を下げるため）
    start_list <- list(
      base_start,
      c(log(0.8), log(K0 * 1.5), log(0.02), log(B1_approx), log(max(0.1, sigma0))),
      c(log(0.3), log(K0 * 0.7), log(0.1), log(B1_approx), log(max(0.1, sigma0)))
    )
  }

  best <- list(nll = Inf, par = NA, converged = FALSE)

  for (st in start_list) {
    fit <- tryCatch(
      optim(st, nll_fun, method = "BFGS", control = control),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < best$nll) {
      best$nll <- fit$value
      best$par <- fit$par
      best$converged <- (fit$convergence == 0)
    }
  }

  if (!is.finite(best$nll)) {
    return(list(
      success = FALSE,
      message = "optim failed",
      par = NA,
      B_hat = rep(NA_real_, TTYear),
      nll = NA_real_,
      converged = FALSE
    ))
  }

  est <- c(
    r = exp(best$par[1]),
    K = exp(best$par[2]),
    q = exp(best$par[3]),
    B1 = exp(best$par[4]),
    sigma_I = exp(best$par[5])
  )

  B_hat <- predict_B_schaefer(C_obs, est["r"], est["K"], est["B1"])
  if (any(!is.finite(B_hat))) {  # 変更点: 推定後B_hatの健全性チェック
    return(list(
      success = FALSE,
      message = "B_hat is not finite under estimated parameters",
      par = est,
      B_hat = B_hat,
      nll = best$nll,
      converged = best$converged
    ))
  }

  list(
    success = TRUE,
    par = est,
    B_hat = B_hat,
    nll = best$nll,
    converged = best$converged
  )
}

# ----------------------------------------------------------
# C++（Rcpp）版：最尤の計算を高速化
# ----------------------------------------------------------
# Rcppが作成する関数の格納先（空環境だと sourceCpp 内で <-
# が見えずに失敗するため、親は globalenv にする）
.cpp_env <- new.env(parent = globalenv())

fit_schaefer_ml_cpp <- function(C_obs, CPUE_obs, start_list = NULL, control = list(maxit = 1000),
                                fallback_to_r = TRUE) {
  # Rcppがない場合はR版にフォールバック
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    message("Rcpp is not available. Use R version.")
    if (fallback_to_r) return(fit_schaefer_ml(C_obs, CPUE_obs, start_list = start_list, control = control))
    stop("Rcpp is not available.")
  }

  TTYear <- length(C_obs)
  if (length(CPUE_obs) != TTYear) stop("length mismatch: C_obs and CPUE_obs")

  # CPUEの正値チェック（対数正規の前提）
  idx <- which(is.finite(CPUE_obs) & CPUE_obs > 0)
  if (length(idx) < 5) {
    return(list(
      success = FALSE,
      message = "CPUE_obs has too few positive values",
      par = NA,
      B_hat = rep(NA_real_, TTYear),
      nll = NA_real_,
      converged = FALSE
    ))
  }

  # C++関数を一度だけコンパイル
  if (!exists("nll_schaefer_cpp", envir = .cpp_env, inherits = FALSE)) {
    # cppFunctionは1関数のみ対応のため、sourceCppでまとめて定義
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      using namespace Rcpp;

      // [[Rcpp::export]]
      NumericVector predict_B_schaefer_cpp(NumericVector C_obs, double r, double K, double B1) {
        int TT = C_obs.size();
        NumericVector B(TT);
        B[0] = B1;
        if (TT >= 2) {
          for (int t = 0; t < TT - 1; ++t) {
            double B_next = B[t] + r * B[t] * (1.0 - B[t] / K) - C_obs[t];
            if (!R_finite(B_next)) return rep(NA_REAL, TT);
            if (B_next <= 0.0) return rep(NA_REAL, TT);
            B[t + 1] = B_next;
          }
        }
        return B;
      }

      // [[Rcpp::export]]
      double nll_schaefer_cpp(NumericVector par, NumericVector C_obs, NumericVector CPUE_obs,
                              IntegerVector idx) {
        if (par.size() != 5) return 1e12; // 変更点: par長チェック
        double r = std::exp(par[0]);
        double K = std::exp(par[1]);
        double q = std::exp(par[2]);
        double B1 = std::exp(par[3]);
        double sigma_I = std::exp(par[4]);

        int TT = C_obs.size();
        NumericVector B(TT);
        B[0] = B1;
        for (int t = 0; t < TT - 1; ++t) {
          double B_next = B[t] + r * B[t] * (1.0 - B[t] / K) - C_obs[t];
          if (!R_finite(B_next)) return 1e12;
          if (B_next <= 0.0) return 1e12;
          B[t + 1] = B_next;
        }

        double nll = 0.0;
        double log2pi = std::log(2.0 * 3.14159265358979323846);
        for (int i = 0; i < idx.size(); ++i) {
          int t = idx[i] - 1; // 1-based -> 0-based
          if (t < 0 || t >= TT) continue;
          double B_t = B[t];
          if (!R_finite(B_t) || B_t <= 0) return 1e12;
          double mu = std::log(q * B_t) - 0.5 * sigma_I * sigma_I;
          double x = std::log(CPUE_obs[t]);
          double z = (x - mu) / sigma_I;
          nll += 0.5 * (log2pi + 2.0 * std::log(sigma_I) + z * z);
        }
        return nll;
      }
    ', env = .cpp_env, verbose = FALSE)
  }

  # 負の対数尤度（C++）呼び出し
  idx_cpp <- as.integer(idx)
  nll_fun <- function(par) {
    if (length(par) != 5L) return(1e12)  # 変更点: 安全チェック
    val <- .cpp_env$nll_schaefer_cpp(par, C_obs, CPUE_obs, idx_cpp)
    if (!is.finite(val)) 1e12 else val
  }

  # 初期値（R版と同じ）
  q0 <- 0.05
  first <- idx[1]
  B1_approx <- max(CPUE_obs[first] / q0, 1e-6)
  K0 <- max(B1_approx * 2, max(C_obs, na.rm = TRUE) * 5, 1)
  r0 <- 0.5
  sigma0 <- max(0.05, sd(log(CPUE_obs[idx]), na.rm = TRUE))
  if (!is.finite(sigma0)) sigma0 <- 0.2

  base_start <- c(log(r0), log(K0), log(q0), log(B1_approx), log(sigma0))
  if (is.null(start_list)) {
    start_list <- list(
      base_start,
      c(log(0.8), log(K0 * 1.5), log(0.02), log(B1_approx), log(max(0.1, sigma0))),
      c(log(0.3), log(K0 * 0.7), log(0.1), log(B1_approx), log(max(0.1, sigma0)))
    )
  }

  best <- list(nll = Inf, par = NA, converged = FALSE)
  for (st in start_list) {
    if (length(st) != 5L) next
    
    fit <- tryCatch(
      optim(st, nll_fun, method = "BFGS", control = control),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < best$nll) {
      best$nll <- fit$value
      best$par <- fit$par
      best$converged <- (fit$convergence == 0)
    }
  }

  if (!is.finite(best$nll)) {
    return(list(
      success = FALSE,
      message = "optim failed",
      par = NA,
      B_hat = rep(NA_real_, TTYear),
      nll = NA_real_,
      converged = FALSE
    ))
  }

  est <- c(
    r = exp(best$par[1]),
    K = exp(best$par[2]),
    q = exp(best$par[3]),
    B1 = exp(best$par[4]),
    sigma_I = exp(best$par[5])
  )

  B_hat <- .cpp_env$predict_B_schaefer_cpp(C_obs, est["r"], est["K"], est["B1"])
  if (any(!is.finite(B_hat))) {  # 変更点: 推定後B_hatの健全性チェック
    return(list(
      success = FALSE,
      message = "B_hat is not finite under estimated parameters",
      par = est,
      B_hat = B_hat,
      nll = best$nll,
      converged = best$converged
    ))
  }

  list(
    success = TRUE,
    par = est,
    B_hat = B_hat,
    nll = best$nll,
    converged = best$converged
  )
}
