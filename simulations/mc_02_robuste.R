# =============================================================================
# Monte Carlo — Module 2 : hétéroscédasticité, couverture et robustesse
# Sous un DGP hétéroscédastique connu, on montre que :
#  - l'OLS reste sans biais ;
#  - les IC classiques (s^2 (X'X)^{-1}) SOUS-COUVRENT ;
#  - les IC robustes HC3 (et HC0) RESTAURENT la couverture ~0.95 ;
#  - le test t classique sur un coef nul rejette trop souvent (> 5 %),
#    HC3 ramène le taux de rejet au niveau nominal.
#
# DGP : y = b0 + b1 x1 + b2 x2 + eps, eps_i ~ N(0, sigma_i^2),
#       sigma_i = exp(gamma * x2_i)  (variance pilotée par x2).
#       b2 = 0 : le coef nul suivi est CELUI dont le régresseur pilote la
#       variance -> distorsion maximale du test classique sur b2.
# =============================================================================

for (f in c("00_linalg", "01_ols", "02_gls_robust"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- DGP (X fixe) -----------------------------------------------------------
n <- 40L
beta_true <- c(b0 = 1, b1 = 2, b2 = 0)
gamma <- 0.8
x1 <- rnorm(n); x2 <- rnorm(n)
sigma_i <- exp(gamma * x2)                 # variance pilotée par x2 (coef nul)
dat <- data.frame(x1 = x1, x2 = x2)
Xmat <- model.matrix(~ x1 + x2, dat)
df_res <- n - 3L
t95 <- qt(0.975, df_res)

R <- 5000L
p <- 3L
beta_hat <- matrix(NA_real_, R, p)
cover <- list(classique = matrix(FALSE, R, p),
              HC0 = matrix(FALSE, R, p),
              HC3 = matrix(FALSE, R, p))
reject_b2 <- c(classique = 0, HC0 = 0, HC3 = 0)   # rejet de H0: b2 = 0 à 5 %

for (r in seq_len(R)) {
  dat$y <- as.numeric(Xmat %*% beta_true + rnorm(n, sd = sigma_i))
  fit <- ols_fit(y ~ x1 + x2, dat)
  b <- fit$coefficients
  beta_hat[r, ] <- b
  se_list <- list(classique = sqrt(diag(fit$vcov)),
                  HC0 = sqrt(diag(vcov_hc(fit, "HC0"))),
                  HC3 = sqrt(diag(vcov_hc(fit, "HC3"))))
  for (m in names(se_list)) {
    se <- se_list[[m]]
    cover[[m]][r, ] <- (b - t95 * se <= beta_true) & (beta_true <= b + t95 * se)
    # test H0: b2 = 0
    if (abs(b[3] / se[3]) > t95) reject_b2[m] <- reject_b2[m] + 1
  }
}

# ---- Résultats --------------------------------------------------------------
cat("=== Non-biais : moyenne(beta_hat) - beta_vrai ===\n")
print(round(colMeans(beta_hat) - beta_true, 4))

cov_tab <- data.frame(
  coef = names(beta_true),
  classique = round(colMeans(cover$classique), 3),
  HC0 = round(colMeans(cover$HC0), 3),
  HC3 = round(colMeans(cover$HC3), 3))
cat("\n=== Couverture empirique des IC à 95 % (R =", R, ") ===\n")
print(cov_tab, row.names = FALSE)

cat("\n=== Taux de rejet de H0: b2 = 0 (niveau nominal 5 %) ===\n")
print(round(reject_b2 / R, 3))
cat("Classique -> gonflé (>5 %) ; HC3 -> ~5 % (correct).\n")

# ---- Graphique couverture ---------------------------------------------------
long <- do.call(rbind, lapply(c("classique", "HC0", "HC3"), function(m)
  data.frame(coef = names(beta_true), couv = colMeans(cover[[m]]), methode = m)))
long$methode <- factor(long$methode, levels = c("classique", "HC0", "HC3"))

gg <- ggplot(long, aes(coef, couv, fill = methode)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.80, 1)) +
  labs(title = "Couverture des IC à 95 % sous hétéroscédasticité",
       subtitle = "Classique : sous-couvre (surtout b2, dont le régresseur pilote la variance) ; HC3 restaure ~0.95",
       x = NULL, y = "couverture empirique", fill = "Variance") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_02_couverture.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "mc_02_couverture.png"), "\n")

# ---- Interprétation ---------------------------------------------------------
# - beta_hat sans biais : l'hétéroscédasticité n'affecte pas le point estimate.
# - IC classiques : la couverture tombe sous 0.95 (Prop. 2.1) ; l'effet est le
#   plus marqué sur b2 (le coef dont le régresseur x2 pilote la variance).
# - HC3 (et HC0) : couverture ramenée vers 0.95, et taux de rejet sous H0
#   ramené vers 5 %. HC3 est le plus fiable en petit échantillon.
