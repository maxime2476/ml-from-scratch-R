# =============================================================================
# Monte Carlo — Module 18 : GMM
#  (A) Efficacité : sous hétéroscédasticité, la GMM efficace (2 étapes) a une
#      variance PLUS FAIBLE que le 2SLS pour le coefficient endogène.
#  (B) Test de suridentification J : taille ~5% sous instruments valides ;
#      puissance quand un instrument est invalide (corrélé à l'erreur).
# =============================================================================

for (f in c("00_linalg", "05_iv_2sls", "18_gmm", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

b_end <- -1.5
gen <- function(n, invalid = 0) {
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n); x1 <- rnorm(n)
  xend <- 0.9 * z1 + 0.7 * z2 + w + rnorm(n)
  u <- 2 * w + rnorm(n) * exp(0.5 * z1)             # hétéroscédastique
  # instrument éventuellement INVALIDE : z2 corrélé à u si invalid>0
  z2 <- z2 + invalid * u
  y <- 1 + 2 * x1 + b_end * xend + u
  list(y = y, X = cbind(1, x1, xend), Z = cbind(1, x1, z1, z2))
}

# =============================================================================
# (A) Efficacité 2SLS vs GMM efficace (coef endogène)
# =============================================================================
n <- 400; R <- 2000
b2sls <- bgmm <- numeric(R); cov2 <- covg <- logical(R)
for (r in seq_len(R)) {
  d <- gen(n)
  g1 <- gmm_linear(d$y, d$X, d$Z, twostep = FALSE)
  g2 <- gmm_linear(d$y, d$X, d$Z, twostep = TRUE)
  b2sls[r] <- g1$coefficients[3]; bgmm[r] <- g2$coefficients[3]
  cov2[r] <- abs(g1$coefficients[3] - b_end) <= 1.96 * g1$se[3]
  covg[r] <- abs(g2$coefficients[3] - b_end) <= 1.96 * g2$se[3]
}
s2 <- mc_summary(b2sls, b_end); sg <- mc_summary(bgmm, b_end)
cat("=== (A) Efficacité (coef endogène, vrai =", b_end, ", n =", n, ") ===\n")
cat(sprintf("2SLS         : RMSE %.4f (%.4f)  sd empirique %.4f  couverture %.3f\n",
            s2$rmse, s2$rmse_se, s2$empirical_se, mean(cov2)))
cat(sprintf("GMM efficace : RMSE %.4f (%.4f)  sd empirique %.4f  couverture %.3f\n",
            sg$rmse, sg$rmse_se, sg$empirical_se, mean(covg)))
cat(sprintf("Gain de variance GMM/2SLS : %.1f %%\n", 100 * (1 - sg$empirical_se^2 / s2$empirical_se^2)))

dfA <- rbind(data.frame(b = b2sls, est = "2SLS"), data.frame(b = bgmm, est = "GMM efficace"))
ggA <- ggplot(dfA, aes(b, fill = est)) + geom_density(alpha = 0.5) +
  geom_vline(xintercept = b_end, linetype = "dashed") +
  labs(title = "GMM efficace vs 2SLS sous hétéroscédasticité",
       subtitle = "La GMM efficace (pondération S^-1) concentre davantage : variance plus faible",
       x = expression(hat(b)[endogene]), y = "densité", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_18_efficacite.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Test J : taille et puissance
# =============================================================================
cat("\n=== (B) Test de suridentification J (niveau 5 %) ===\n")
for (inv in c(0, 0.3, 0.6)) {
  rej <- mean(replicate(1500, {
    d <- gen(300, invalid = inv)
    gmm_linear(d$y, d$X, d$Z, twostep = TRUE)$J_pvalue < 0.05
  }))
  rm <- reject_mc(rep(c(TRUE, FALSE), c(round(rej * 1500), 1500 - round(rej * 1500))))
  lbl <- if (inv == 0) "instruments valides (taille)" else sprintf("z2 invalide (%.1f) — puissance", inv)
  cat(sprintf("  %-32s : rejet %.3f ± %.3f\n", lbl, rm$rate, rm$se))
}
cat("Sous instruments valides : ~5 %. Instrument invalide : le J le détecte (puissance).\n")
cat("\nGraphique -> ", file.path(out_dir, "mc_18_efficacite.png"), "\n")
