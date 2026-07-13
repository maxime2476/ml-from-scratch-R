# =============================================================================
# Monte Carlo — Module 20 : régression quantile
#  (A) Sous hétéroscédasticité (modèle location-scale), les pentes des quantiles
#      DIFFÈRENT ("éventail") : la régression quantile les retrouve, l'OLS ne voit
#      qu'une pente moyenne constante.
#  (B) Robustesse : la régression médiane (LAD) résiste aux valeurs aberrantes
#      (RMSE stable), l'OLS est cassé.
# =============================================================================

for (f in c("00_linalg", "01_ols", "20_quantile", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# DGP location-scale : y = 1 + 2 x + (1 + gamma*x) * eps, eps ~ N(0,1), x>0.
# Pente vraie du tau-quantile : 2 + gamma * qnorm(tau).
gamma <- 0.8
true_slope <- function(tau) 2 + gamma * qnorm(tau)

# =============================================================================
# (A) Éventail des pentes de quantiles
# =============================================================================
n <- 500; R <- 500
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
est <- matrix(NA_real_, R, length(taus), dimnames = list(NULL, as.character(taus)))
ols_slope <- numeric(R)
for (r in seq_len(R)) {
  x <- runif(n, 0.2, 2); y <- 1 + 2 * x + (1 + gamma * x) * rnorm(n)
  d <- data.frame(x = x, y = y)
  for (j in seq_along(taus)) est[r, j] <- qreg_fit(y ~ x, d, tau = taus[j])$coefficients["x"]
  ols_slope[r] <- ols_fit_slope <- coef(lm(y ~ x, d))["x"]
}
cat("=== (A) Pentes des quantiles (vrai = 2 + 0.8*qnorm(tau)) ===\n")
tab <- data.frame(tau = taus, vrai = round(true_slope(taus), 3),
                  estime = round(colMeans(est), 3),
                  mc_se = round(apply(est, 2, mc_se), 4))
print(tab, row.names = FALSE)
cat(sprintf("OLS : pente moyenne = %.3f (constante — ne capte pas l'éventail)\n", mean(ols_slope)))

dfA <- data.frame(tau = taus, estime = colMeans(est),
                  lo = colMeans(est) - 1.96 * apply(est, 2, mc_se),
                  hi = colMeans(est) + 1.96 * apply(est, 2, mc_se),
                  vrai = true_slope(taus))
ggA <- ggplot(dfA, aes(tau)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "#2166ac") +
  geom_line(aes(y = estime, colour = "quantile estimé"), linewidth = 1) +
  geom_point(aes(y = estime, colour = "quantile estimé"), size = 2) +
  geom_line(aes(y = vrai, colour = "quantile vrai"), linetype = "dashed") +
  geom_hline(aes(yintercept = mean(ols_slope), colour = "OLS (moyenne)"), linetype = "dotted") +
  labs(title = "Régression quantile : éventail des pentes sous hétéroscédasticité",
       subtitle = "La pente varie avec tau (location-scale) ; l'OLS ne voit qu'une pente moyenne",
       x = expression(tau), y = "pente de x", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_20_eventail.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Robustesse : LAD vs OLS sous contamination
# =============================================================================
cat("\n=== (B) Robustesse : RMSE de la pente (vrai = 2) selon %% d'outliers ===\n")
Rb <- 800; nb <- 200
for (frac in c(0, 0.05, 0.15)) {
  bo <- bl <- numeric(Rb)
  for (r in seq_len(Rb)) {
    x <- rnorm(nb); y <- 1 + 2 * x + rnorm(nb)
    k <- floor(frac * nb); if (k > 0) y[seq_len(k)] <- y[seq_len(k)] + 30
    d <- data.frame(x = x, y = y)
    bo[r] <- coef(lm(y ~ x, d))["x"]; bl[r] <- qreg_fit(y ~ x, d, 0.5)$coefficients["x"]
  }
  cat(sprintf("  %4.0f%% outliers : RMSE(OLS) = %.3f   RMSE(LAD) = %.3f\n",
              100 * frac, mc_summary(bo, 2)$rmse, mc_summary(bl, 2)$rmse))
}
cat("La médiane (LAD) garde une RMSE faible ; l'OLS explose avec les outliers.\n")
cat("\nGraphique -> ", file.path(out_dir, "mc_20_eventail.png"), "\n")
