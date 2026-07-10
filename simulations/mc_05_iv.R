# =============================================================================
# Monte Carlo — Module 5 : variables instrumentales
#  (A) Endogénéité : l'OLS est biaisé et ses IC sous-couvrent ; le 2SLS est
#      consistant et ses IC couvrent ~0.95.
#  (B) Instruments faibles : quand la force des instruments (F de 1ʳᵉ étape)
#      diminue, le 2SLS se biaise vers l'OLS et sa couverture s'effondre.
#
# DGP : x_end = pi*(z1+z2) + rho_w * w + nu ; u = rho_w * w + e (confondeur w) ;
#       y = b0 + b1*x1 + b2*x_end + u.  Cov(x_end, u) != 0 via w -> endogénéité.
# =============================================================================

for (f in c("00_linalg", "01_ols", "05_iv_2sls"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

b0 <- 1; b1 <- 2; b2 <- -1.5
t95 <- function(df) qt(0.975, df)

gen <- function(n, pi_str) {
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n); x1 <- rnorm(n)
  x_end <- pi_str * (z1 + z2) + 2.0 * w + rnorm(n)
  u <- 2.0 * w + rnorm(n)                        # forte corrélation avec x_end via w
  y <- b0 + b1 * x1 + b2 * x_end + u
  list(y = y,
       X = cbind(1, x1, x_end),
       Z = cbind(1, x1, z1, z2))
}

# =============================================================================
# (A) Endogénéité : OLS vs 2SLS (instruments forts)
# =============================================================================
n <- 300; R <- 3000
res <- matrix(NA_real_, R, 4,
              dimnames = list(NULL, c("ols_b2", "tsls_b2", "ols_cov", "tsls_cov")))
for (r in seq_len(R)) {
  d <- gen(n, pi_str = 0.7)
  ols <- solve_ls_qr(d$X, d$y)
  se_ols <- sqrt(sum(ols$residuals^2) / (n - 3) * diag(chol2inv(chol(crossprod(d$X)))))
  b_ols <- ols$coefficients[3]
  fit <- tsls_fit(d$y, d$X, d$Z)
  b_ts <- fit$coefficients[3]; se_ts <- fit$se[3]
  res[r, ] <- c(b_ols, b_ts,
                abs(b_ols - b2) <= t95(n - 3) * se_ols[3],
                abs(b_ts  - b2) <= t95(fit$df.residual) * se_ts)
}
cat("=== (A) Endogénéité (instruments forts, n =", n, ") ===\n")
cat(sprintf("Biais OLS  (b2) : %+.4f   couverture IC : %.3f\n",
            mean(res[, "ols_b2"]) - b2, mean(res[, "ols_cov"])))
cat(sprintf("Biais 2SLS (b2) : %+.4f   couverture IC : %.3f\n",
            mean(res[, "tsls_b2"]) - b2, mean(res[, "tsls_cov"])))
cat("OLS : biaisé + sous-couverture ; 2SLS : ~sans biais + couverture ~0.95.\n")

dens <- rbind(data.frame(b2 = res[, "ols_b2"], est = "OLS"),
              data.frame(b2 = res[, "tsls_b2"], est = "2SLS"))
ggA <- ggplot(dens, aes(b2, fill = est)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = b2, linetype = "dashed") +
  labs(title = "Endogénéité : loi d'échantillonnage de l'estimateur de b2",
       subtitle = "OLS centré loin de la vraie valeur (pointillé) ; 2SLS centré dessus",
       x = expression(hat(b)[2]), y = "densité", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_05_endogeneite.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Instruments faibles : force pi décroissante
# =============================================================================
pis <- c(0.5, 0.25, 0.12, 0.06, 0.03)
Rb <- 2000; nb <- 300
tabB <- data.frame(pi = pis, F_moyen = NA, biais_2sls = NA, couv_2sls = NA, biais_ols = NA)
for (i in seq_along(pis)) {
  pu <- pis[i]
  bs <- numeric(Rb); cov <- 0; Fs <- numeric(Rb); bo <- numeric(Rb)
  for (r in seq_len(Rb)) {
    d <- gen(nb, pi_str = pu)
    fit <- tsls_fit(d$y, d$X, d$Z)
    bs[r] <- fit$coefficients[3]
    cov <- cov + (abs(fit$coefficients[3] - b2) <= t95(fit$df.residual) * fit$se[3])
    Fs[r] <- first_stage_F(d$X[, 3], d$Z, excluded = c(3, 4))$F
    bo[r] <- solve_ls_qr(d$X, d$y)$coefficients[3]
  }
  tabB[i, -1] <- c(mean(Fs), mean(bs) - b2, cov / Rb, mean(bo) - b2)
}
cat("\n=== (B) Instruments faibles (n =", nb, ") ===\n")
print(round(tabB, 3), row.names = FALSE)
cat("Quand F chute (< 10), le 2SLS se biaise vers l'OLS et la couverture s'effondre.\n")

ggB <- ggplot(tabB, aes(F_moyen)) +
  geom_line(aes(y = biais_2sls, colour = "biais 2SLS"), linewidth = 1) +
  geom_point(aes(y = biais_2sls, colour = "biais 2SLS")) +
  geom_hline(aes(yintercept = mean(tabB$biais_ols), colour = "biais OLS (réf.)"), linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_x_log10() +
  labs(title = "Instruments faibles : biais du 2SLS vs force des instruments",
       subtitle = "Quand le F de 1ʳᵉ étape chute, le 2SLS dérive vers le biais de l'OLS",
       x = "F de première étape (moyen, échelle log)", y = "biais de b2", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_05_instruments_faibles.png"), ggB, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_05_endogeneite.png, mc_05_instruments_faibles.png\n", sep = "")
