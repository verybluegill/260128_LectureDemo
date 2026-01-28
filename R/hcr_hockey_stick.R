# -*- coding: utf-8 -*-
# (saved as UTF-8, no BOM)

# ホッケースティックHCR（資源量ベース）
hcr_hockey_stick <- function(B_hat, beta = 0.8, Bban, Blim) {
  if (Blim <= Bban) stop("Blim must be > Bban")

  B_hat <- pmax(B_hat, 0)
  h <- rep(NA_real_, length(B_hat))

  h[B_hat <= Bban] <- 0
  h[B_hat >= Blim] <- beta

  mid <- B_hat > Bban & B_hat < Blim
  h[mid] <- beta * (B_hat[mid] - Bban) / (Blim - Bban)

  C_next <- h * B_hat

  list(h = h, C_next = C_next)
}

# 形状の図示用（Bとhの対応）
hcr_hockey_shape <- function(B_seq, beta = 0.8, Bban, Blim) {
  out <- hcr_hockey_stick(B_seq, beta = beta, Bban = Bban, Blim = Blim)
  data.frame(B = B_seq, h = out$h)
}