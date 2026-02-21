# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)
# 本当はHockey-stickのHCRはFmsyを参照すべき!
# ホッケースティックHCR（資源量ベース）
hcr_hockey_stick <- function(B_hat, beta = 0.8, Bban, Blim, Fmsy = 1) {
  if (Blim <= Bban) stop("Blim must be > Bban")
  
  if (any(!is.finite(B_hat))) stop("B_hat contains non-finite (Inf/NaN)")  
  if (!is.finite(beta) || beta < 0) stop("beta must be non-negative")      
  if (!is.finite(Bban) || !is.finite(Blim)) stop("Bban/Blim must be finite") 
  if (!is.finite(Fmsy) || Fmsy < 0) stop("Fmsy must be non-negative")
  
  B_hat <- pmax(B_hat, 0)  # 負の推定値は0に丸める
  
  h <- rep(NA_real_, length(B_hat))
  
  h[B_hat <= Bban] <- 0
  h[B_hat >= Blim] <- beta
  
  mid <- B_hat > Bban & B_hat < Blim
  h[mid] <- beta * (B_hat[mid] - Bban) / (Blim - Bban)
  
  F_next <- h * Fmsy
  C_next <- F_next * B_hat
  
  list(h = h, F_next = F_next, C_next = C_next)
}

# 形状の図示用（Bとhの対応）
hcr_hockey_shape <- function(B_seq, beta = 0.8, Bban, Blim) {
  out <- hcr_hockey_stick(B_seq, beta = beta, Bban = Bban, Blim = Blim)
  data.frame(B = B_seq, h = out$h)
}
