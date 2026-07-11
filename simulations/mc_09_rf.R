# =============================================================================
# Monte Carlo — Module 9 : bagging et forêts aléatoires
#  (A) Validation de Var(f̄) = rho*sigma² + (1-rho)/B * sigma² (Prop. 9.1) :
#      estimation de sigma² et rho, puis comparaison à la variance empirique du
#      prédicteur agrégé pour une gamme de B.
#  (B) Effet de mtry : réduire mtry décorréle les arbres (rho baisse), donc
#      abaisse le plancher rho*sigma² — au prix d'un peu de biais/variance.
# =============================================================================

for (f in c("00_linalg", "08_cart", "09_bagging_rf"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# DGP de régression ; point de test fixe x0.
f_true <- function(X) sin(2 * X[, 1]) + X[, 2] - 0.5 * X[, 3] * X[, 4]
gen_data <- function(n, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * 4), n, 4)
  data.frame(x1 = X[, 1], x2 = X[, 2], x3 = X[, 3], x4 = X[, 4],
             y = f_true(X) + rnorm(n, sd = 0.5))
}
x0 <- data.frame(x1 = 0.5, x2 = -0.3, x3 = 0.4, x4 = 0.2)

# Ajuste un arbre sur un bootstrap et prédit en x0.
one_tree_pred <- function(data, mtry) {
  ib <- sample.int(nrow(data), nrow(data), replace = TRUE)
  tr <- cart_fit(y ~ x1 + x2 + x3 + x4, data[ib, ], "anova",
                 max_depth = 20, min_split = 5, min_leaf = 1, mtry = mtry)
  predict_cart(tr, x0)
}

# sigma² = variance d'un arbre (sur les jeux) ; c = Cov entre deux arbres du même
# jeu = Var_D(E[arbre|D]), estimée par la covariance de deux moyennes d'arbres
# DISJOINTES (conditionnellement indépendantes) -> estimateur peu bruité de rho.
est_sigma2_rho <- function(Preds) {
  B <- ncol(Preds); h <- B %/% 2L
  s2 <- mean(apply(Preds, 2, var))                       # variance moyenne d'un arbre
  c_hat <- cov(rowMeans(Preds[, seq_len(h), drop = FALSE]),
               rowMeans(Preds[, (h + 1L):(2L * h), drop = FALSE]))
  list(sigma2 = s2, rho = c_hat / s2)
}

# =============================================================================
# (A) Validation de la formule (9.1) — mtry = 4 (bagging)
# =============================================================================
R <- 250; Bmax <- 60; n <- 150; mtry <- 4
Preds <- matrix(NA_real_, R, Bmax)          # R jeux d'apprentissage x Bmax arbres
for (r in seq_len(R)) {
  d <- gen_data(n, seed = 1000 + r)          # nouveau jeu d'apprentissage
  Preds[r, ] <- vapply(seq_len(Bmax), function(b) one_tree_pred(d, mtry), numeric(1))
}
sr <- est_sigma2_rho(Preds); sigma2 <- sr$sigma2; rho <- sr$rho
Bgrid <- c(1, 2, 5, 10, 20, 40, 60)
tabA <- data.frame(B = Bgrid,
                   var_empirique = sapply(Bgrid, function(B) var(rowMeans(Preds[, seq_len(B), drop = FALSE]))),
                   var_formule = rho * sigma2 + (1 - rho) / Bgrid * sigma2)
cat("=== (A) Validation Var(f̄) = rho*sigma² + (1-rho)/B * sigma² ===\n")
cat(sprintf("sigma² (arbre) = %.4f ; rho = %.3f ; plancher rho*sigma² = %.4f\n",
            sigma2, rho, rho * sigma2))
print(round(tabA, 4), row.names = FALSE)
cat("La variance empirique du prédicteur agrégé suit la formule et tend vers le plancher.\n")

longA <- reshape(tabA, varying = c("var_empirique", "var_formule"), v.names = "variance",
                 timevar = "source", times = c("empirique", "formule (9.1)"), direction = "long")
ggA <- ggplot(longA, aes(B, variance, colour = source)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_hline(yintercept = rho * sigma2, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(Bgrid), y = rho * sigma2, label = "plancher rho*sigma²",
           vjust = -0.6, hjust = 1, size = 3) +
  labs(title = "Bagging : variance du prédicteur agrégé selon B",
       subtitle = "Var(f̄) = rho*sigma² + (1-rho)/B * sigma² -> plancher rho*sigma² (Prop. 9.1)",
       x = "nombre d'arbres B", y = "variance de la prédiction en x0", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_09_variance_B.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Effet de mtry : rho, plancher, et erreur de test
# =============================================================================
mtrys <- c(1, 2, 3, 4)
R2 <- 200; Btwo <- 30
tabB <- data.frame(mtry = mtrys, sigma2 = NA, rho = NA, plancher = NA)
for (im in seq_along(mtrys)) {
  P2 <- matrix(NA_real_, R2, Btwo)
  for (r in seq_len(R2)) {
    d <- gen_data(n, seed = 5000 + r)
    P2[r, ] <- vapply(seq_len(Btwo), function(b) one_tree_pred(d, mtrys[im]), numeric(1))
  }
  sr2 <- est_sigma2_rho(P2)
  tabB[im, -1] <- c(sr2$sigma2, sr2$rho, sr2$rho * sr2$sigma2)
}
cat("\n=== (B) Effet de mtry sur la corrélation et le plancher de variance ===\n")
print(round(tabB, 4), row.names = FALSE)
cat("Réduire mtry -> rho plus faible (arbres décorrélés) MAIS sigma² plus élevé\n",
    "(chaque arbre voit moins de variables) : le plancher rho*sigma² reflète cet\n",
    "arbitrage ; l'erreur de test est optimale à un mtry intermédiaire.\n", sep = "")

# Erreur de test d'une forêt selon mtry (B fixe)
dtrain <- gen_data(400, seed = 42); dtest <- gen_data(3000, seed = 43)
errB <- data.frame(mtry = mtrys, test_mse = NA, oob = NA)
for (im in seq_along(mtrys)) {
  fit <- random_forest_fit(y ~ x1 + x2 + x3 + x4, dtrain, "anova", B = 150,
                           mtry = mtrys[im], seed = 1)
  errB$test_mse[im] <- mean((predict_forest(fit, dtest) - dtest$y)^2)
  errB$oob[im] <- fit$oob_error
}
cat("\n=== Erreur de test / OOB de la forêt selon mtry (B = 150) ===\n")
print(round(errB, 4), row.names = FALSE)

ggB <- ggplot(tabB, aes(mtry)) +
  geom_line(aes(y = rho, colour = "rho (corrélation)"), linewidth = 1) +
  geom_point(aes(y = rho, colour = "rho (corrélation)"), size = 2) +
  geom_line(aes(y = plancher, colour = "plancher rho*sigma²"), linewidth = 1) +
  geom_point(aes(y = plancher, colour = "plancher rho*sigma²"), size = 2) +
  labs(title = "Forêt aléatoire : mtry décorréle les arbres",
       subtitle = "rho croît avec mtry (moins de décorrélation) ; le plancher reflète l'arbitrage avec sigma²",
       x = "mtry (variables candidates par split)", y = "valeur", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_09_mtry.png"), ggB, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_09_variance_B.png, mc_09_mtry.png\n", sep = "")
