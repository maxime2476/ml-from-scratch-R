# =============================================================================
# Monte Carlo — Module 14 : où l'inférence valide casse
#  (A) Inférence post-sélection (OBLIGATOIRE) : sous H0 (tous les coefficients
#      nuls), le t-test OLS post-lasso rejette BIEN AU-DELÀ de 5 % ; l'OLS complet
#      et le sample-splitting restent à ~5 %.
#  (B) Biais de régularisation : les IC centrés sur le ridge SOUS-COUVRENT la
#      vraie valeur ; l'OLS (sans biais) couvre à ~95 %.
# =============================================================================

for (f in c("00_linalg", "01_ols", "04_regularisation", "14_m_estimation"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) Inférence post-sélection sous H0 (aucun signal)
# =============================================================================
n <- 100; p <- 30; R <- 2000; alpha <- 0.05
crit <- qt(1 - alpha / 2, df = n - 1)          # approx ; df ajusté par sous-modèle

t_reject_naive <- 0; n_selected <- 0
t_reject_full  <- 0; n_full <- 0
t_reject_split <- 0; n_split <- 0

lasso_select <- function(Xs, yc) {
  lam <- 0.5 * max(abs(crossprod(Xs, yc)))       # lambda adaptatif -> quelques variables
  b <- lasso_fit(Xs, yc, lambda = lam, standardize = FALSE, intercept = FALSE)$beta
  which(b != 0)
}
ols_pvals <- function(Xsub, y) {
  d <- as.data.frame(Xsub); names(d) <- paste0("v", seq_len(ncol(Xsub))); d$y <- y
  fit <- ols_fit(as.formula(paste("y ~", paste(names(d)[-ncol(d)], collapse = "+"))), d)
  sm <- ols_summary(fit)$coefficients
  sm$p_value[rownames(sm) != "(Intercept)"]
}

for (r in seq_len(R)) {
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)                                  # H0 : y indépendant de X
  Xs <- scale(X); yc <- y - mean(y)

  # --- Post-lasso naïf : sélection ET inférence sur TOUTES les données ---
  S <- lasso_select(Xs, yc)
  if (length(S) > 0 && length(S) < n - 2) {
    pv <- ols_pvals(X[, S, drop = FALSE], y)
    t_reject_naive <- t_reject_naive + sum(pv < alpha); n_selected <- n_selected + length(S)
  }

  # --- OLS complet (sans sélection) : inférence valide ---
  pv_full <- ols_pvals(X, y)
  t_reject_full <- t_reject_full + sum(pv_full < alpha); n_full <- n_full + p

  # --- Sample-splitting : sélection sur moitié A, inférence sur moitié B ---
  iA <- sample.int(n, n / 2); iB <- setdiff(seq_len(n), iA)
  SsA <- lasso_select(scale(X[iA, ]), y[iA] - mean(y[iA]))
  if (length(SsA) > 0 && length(SsA) < length(iB) - 2) {
    pv_split <- ols_pvals(X[iB, SsA, drop = FALSE], y[iB])
    t_reject_split <- t_reject_split + sum(pv_split < alpha); n_split <- n_split + length(SsA)
  }
}

cat("=== (A) Inférence post-sélection sous H0 (niveau nominal 5 %) ===\n")
cat(sprintf("Post-lasso NAÏF  : taux de rejet = %.1f %% (sur %d coefs sélectionnés)\n",
            100 * t_reject_naive / n_selected, n_selected))
cat(sprintf("OLS complet      : taux de rejet = %.1f %% (valide)\n",
            100 * t_reject_full / n_full))
cat(sprintf("Sample-splitting : taux de rejet = %.1f %% (valide)\n",
            100 * t_reject_split / n_split))
cat("=> réutiliser les données pour sélectionner PUIS tester casse l'inférence.\n")

dfA <- data.frame(
  methode = factor(c("post-lasso naïf", "OLS complet", "sample-splitting"),
                   levels = c("post-lasso naïf", "OLS complet", "sample-splitting")),
  rejet = c(t_reject_naive / n_selected, t_reject_full / n_full, t_reject_split / n_split))
ggA <- ggplot(dfA, aes(methode, rejet, fill = methode)) +
  geom_col(width = 0.6) + geom_hline(yintercept = 0.05, linetype = "dashed") +
  annotate("text", x = 0.6, y = 0.06, label = "5 % nominal", hjust = 0, size = 3) +
  labs(title = "Inférence post-sélection cassée (sous H0 : aucun signal)",
       subtitle = "Le t-test post-lasso naïf sur-rejette massivement ; OLS complet et split restent à ~5 %",
       x = NULL, y = "taux de rejet à 5 %", fill = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "mc_14_post_selection.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Biais de régularisation : couverture des IC (ridge vs OLS)
# =============================================================================
n2 <- 80; p2 <- 6; R2 <- 3000
beta_true <- c(2, -1.5, 1, 0.5, -0.8, 0.3); sigma <- 1
lam <- 6
X2 <- matrix(rnorm(n2 * p2), n2, p2)          # design fixe
cov_ols <- cov_ridge <- rep(0, p2)
for (r in seq_len(R2)) {
  y <- as.numeric(X2 %*% beta_true + rnorm(n2, sd = sigma))
  # OLS (sans biais)
  ols <- solve_ls_qr(X2, y); e <- ols$residuals; s2 <- sum(e^2) / (n2 - p2)
  Rinv <- backsolve(ols$R, diag(p2)); se_ols <- sqrt(s2 * rowSums(Rinv^2))
  cov_ols <- cov_ols + (abs(ols$coefficients - beta_true) <= 1.96 * se_ols)
  # Ridge (biaisé) : IC centré sur beta_ridge, variance fréquentiste plug-in
  A <- solve(crossprod(X2) + lam * diag(p2))
  b_ridge <- as.numeric(A %*% crossprod(X2, y))
  V_ridge <- s2 * A %*% crossprod(X2) %*% A          # var fréquentiste du ridge
  se_ridge <- sqrt(diag(V_ridge))
  cov_ridge <- cov_ridge + (abs(b_ridge - beta_true) <= 1.96 * se_ridge)
}
tabB <- data.frame(coef = paste0("b", seq_len(p2)),
                   couv_OLS = round(cov_ols / R2, 3),
                   couv_ridge = round(cov_ridge / R2, 3))
cat("\n=== (B) Couverture des IC à 95 % : OLS (sans biais) vs ridge (biaisé) ===\n")
print(tabB, row.names = FALSE)
cat("Le ridge SOUS-COUVRE : son IC est centré sur un estimateur biaisé (rétréci vers 0).\n")

longB <- rbind(data.frame(coef = tabB$coef, couv = tabB$couv_OLS, est = "OLS"),
               data.frame(coef = tabB$coef, couv = tabB$couv_ridge, est = "ridge"))
ggB <- ggplot(longB, aes(coef, couv, fill = est)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Biais de régularisation : couverture des IC à 95 %",
       subtitle = "OLS couvre ~0.95 ; les IC du ridge sous-couvrent (centrés sur un estimateur biaisé)",
       x = NULL, y = "couverture empirique", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_14_biais_ridge.png"), ggB, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_14_post_selection.png, mc_14_biais_ridge.png\n", sep = "")
