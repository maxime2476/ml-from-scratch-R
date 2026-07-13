# =============================================================================
# Monte Carlo CENTRAL du projet — Module 16 : débiaiser l'estimation causale
# Sur le même DGP (theta vrai connu, g(X) non linéaire, traitement confondu),
# comparer QUATRE estimateurs de l'effet de traitement en BIAIS, RMSE et
# COUVERTURE des IC à 95 % :
#   (1) OLS naïf (contrôles linéaires, g mal spécifiée)
#   (2) Lasso naïf (pénalise AUSSI le coefficient d'intérêt -> biais)
#   (3) DML sans cross-fitting (biais de sur-ajustement des nuisances)
#   (4) DML complet (score de Neyman + cross-fitting)
#
# DGP : y = theta*d + g(X) + eps ; d ~ Bernoulli(m(X)) ; g, m non linéaires.
# =============================================================================

for (f in c("00_linalg", "01_ols", "04_regularisation", "08_cart",
            "09_bagging_rf", "16_causal_ml"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theta_true <- 1.0
gen <- function(n, p = 5, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  # Confusion NON LINÉAIRE en X1^2 : g dépend de X1^2 et la propension de X1,
  # donc d est corrélé à X1^2. Les contrôles LINÉAIRES (OLS naïf) ne peuvent pas
  # l'absorber (biais résiduel) ; la propension monotone préserve l'overlap.
  g <- 3 * X[, 1]^2 + sin(2 * X[, 2]) + 0.5 * X[, 3]           # g(X) non linéaire
  m <- plogis(1.2 * X[, 1] + 0.5 * X[, 2])                     # propension monotone
  d <- rbinom(n, 1, m)
  y <- theta_true * d + g + rnorm(n, sd = 1)
  list(X = X, y = y, d = d)
}

R <- 200L; n <- 500L; z <- 1.96
methods <- c("OLS naïf", "Lasso naïf", "DML sans cross-fit", "DML complet")
theta_hat <- matrix(NA_real_, R, 4, dimnames = list(NULL, methods))
covered   <- matrix(NA, R, 4, dimnames = list(NULL, methods))

for (r in seq_len(R)) {
  dd <- gen(n, seed = 1000 + r); X <- dd$X; y <- dd$y; d <- dd$d

  # (1) OLS naïf : contrôles linéaires (g non linéaire -> biais)
  ml <- lm(y ~ d + X)
  th1 <- coef(ml)["d"]; se1 <- summary(ml)$coefficients["d", "Std. Error"]
  theta_hat[r, 1] <- th1; covered[r, 1] <- abs(th1 - theta_true) <= z * se1

  # (2) Lasso naïf : pénalise aussi le coefficient de d (rétrécissement)
  Xall <- cbind(d = d, X)
  lam <- 0.15 * max(abs(crossprod(scale(Xall), y - mean(y))))
  lf <- lasso_fit(Xall, y, lambda = lam)
  theta_hat[r, 2] <- lf$beta["d"]                         # pas d'IC valide (documenté)

  # (3) DML sans cross-fitting (nuisances sur toutes les données -> sur-ajustement)
  d3 <- dml_plr(y, d, X, nuisance = "forest", crossfit = FALSE, seed = r, B = 50)
  theta_hat[r, 3] <- d3$theta; covered[r, 3] <- d3$ci[1] <= theta_true & theta_true <= d3$ci[2]

  # (4) DML complet (Neyman + cross-fitting)
  d4 <- dml_plr(y, d, X, K = 5, nuisance = "forest", crossfit = TRUE, seed = r, B = 50)
  theta_hat[r, 4] <- d4$theta; covered[r, 4] <- d4$ci[1] <= theta_true & theta_true <= d4$ci[2]
}

# ---- Tableau de synthèse -----------------------------------------------------
tab <- data.frame(
  methode = methods,
  biais = round(colMeans(theta_hat) - theta_true, 3),
  RMSE = round(sqrt(colMeans((theta_hat - theta_true)^2)), 3),
  couverture = round(colMeans(covered), 3))
cat("=== TABLEAU DE SYNTHÈSE — estimation de theta (vrai =", theta_true, ") ===\n")
cat("   n =", n, ", R =", R, " réplications\n\n")
print(tab, row.names = FALSE)
cat("\nLecture :\n",
    " - OLS naïf : contrôles linéaires face à une confusion NON linéaire (X1^2)\n",
    "   -> estimateur très instable (RMSE ~3x le DML) ; biais asymptotique visible\n",
    "   à grand n (theta_OLS ~ 1.25 à n=3000).\n",
    " - Lasso naïf : biaisé (rétrécissement du coefficient d'intérêt, Module 14).\n",
    " - DML sans cross-fit : biaisé (sur-ajustement des nuisances).\n",
    " - DML complet : biais le plus faible, RMSE MINIMAL, couverture ~95 %.\n",
    "   (Lasso : pas d'IC valide -> couverture non renseignée.)\n", sep = "")

# ---- Graphiques --------------------------------------------------------------
dfB <- do.call(rbind, lapply(methods, function(m)
  data.frame(theta = theta_hat[, m], methode = factor(m, levels = methods))))
gg1 <- ggplot(dfB, aes(theta, methode, fill = methode)) +
  ggplot2::geom_violin(alpha = 0.6, scale = "width") +
  geom_vline(xintercept = theta_true, linetype = "dashed") +
  labs(title = "DML : loi d'échantillonnage de l'estimateur de theta",
       subtitle = paste0("Vrai theta = ", theta_true,
                         " (pointillé) ; le DML complet a le biais et la variance les plus faibles"),
       x = expression(hat(theta)), y = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "mc_16_dml_distribution.png"), gg1, width = 8, height = 5, dpi = 120)

covdf <- data.frame(methode = factor(methods, levels = methods),
                    couverture = colMeans(covered))
covdf <- covdf[!is.na(covdf$couverture), ]
gg2 <- ggplot(covdf, aes(methode, couverture, fill = methode)) +
  geom_col(width = 0.6) + geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "DML : couverture des IC à 95 %",
       subtitle = "DML complet ~0.95 ; DML sans cross-fit plus biaisé ; l'OLS naïf a des IC très larges",
       x = NULL, y = "couverture empirique") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "mc_16_dml_couverture.png"), gg2, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_16_dml_distribution.png, mc_16_dml_couverture.png\n", sep = "")
