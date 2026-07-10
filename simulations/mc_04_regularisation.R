# =============================================================================
# Monte Carlo — Module 4 : régularisation
#  (A) Ridge : courbes biais² / variance / EQM en fonction de lambda (analytique
#      éq. 4.5, validées par Monte Carlo), et EQM minimale < EQM de l'OLS.
#  (B) Lasso : capacité à retrouver le SUPPORT vrai (sparsité) selon lambda.
# =============================================================================

for (f in c("00_linalg", "04_regularisation"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) Ridge : biais² / variance / EQM analytiques + validation Monte Carlo
# =============================================================================
n <- 60; p <- 6
rho <- 0.9                                             # forte corrélation -> colinéarité
Sig <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
X <- .standardize(matrix(rnorm(n * p), n, p) %*% chol(Sig))$Xs
beta_true <- c(2, -1.5, 1, 0.5, 0, -0.8)
sigma2 <- 1.5
lambdas <- 10^seq(-2, 3, length.out = 40)

# Courbes analytiques (éq. 4.4-4.5)
ana <- t(sapply(lambdas, function(l) {
  bv <- ridge_bias_var(X, beta_true, sigma2, l)
  c(biais2 = bv$bias2, variance = bv$variance, eqm = bv$mse)
}))
ana <- data.frame(lambda = lambdas, ana)
lam_opt <- lambdas[which.min(ana$eqm)]
eqm_ols <- ana$eqm[which.min(abs(lambdas - min(lambdas)))]

cat("=== (A) Ridge biais-variance ===\n")
cat(sprintf("lambda* (EQM min)      = %.3f\n", lam_opt))
cat(sprintf("EQM(lambda*)           = %.4f\n", min(ana$eqm)))
cat(sprintf("EQM(lambda ~ 0, ~OLS)  = %.4f\n", eqm_ols))
cat(sprintf("gain relatif EQM       = %.1f %%\n", 100 * (1 - min(ana$eqm) / eqm_ols)))

# Validation Monte Carlo sur quelques lambda
check_l <- c(0.1, lam_opt, 50)
R <- 5000
cat("\nValidation Monte Carlo (EQM analytique vs empirique) :\n")
for (l in check_l) {
  B <- matrix(0, R, p)
  for (r in seq_len(R)) {
    y <- as.numeric(X %*% beta_true + rnorm(n, sd = sqrt(sigma2)))
    B[r, ] <- ridge_fit(X, y, l, standardize = FALSE, intercept = FALSE)$beta
  }
  mse_mc <- sum((colMeans(B) - beta_true)^2) + sum(apply(B, 2, var))
  mse_an <- ridge_bias_var(X, beta_true, sigma2, l)$mse
  cat(sprintf("  lambda=%7.3f : EQM analytique=%.4f  MC=%.4f\n", l, mse_an, mse_mc))
}

longA <- reshape(ana, varying = c("biais2", "variance", "eqm"), v.names = "valeur",
                 timevar = "composante", times = c("biais²", "variance", "EQM"),
                 direction = "long")
ggA <- ggplot(longA, aes(lambda, valeur, colour = composante)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = lam_opt, linetype = "dashed", colour = "grey40") +
  scale_x_log10() +
  annotate("text", x = lam_opt, y = max(ana$eqm), label = "lambda*", hjust = -0.1, size = 3.5) +
  labs(title = "Ridge : décomposition biais²-variance de l'EQM",
       subtitle = "Le biais croît, la variance décroît ; l'EQM est minimale en un lambda* > 0",
       x = expression(lambda), y = "valeur", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_04_ridge_bv.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Lasso : récupération du support vrai
# =============================================================================
nB <- 100; pB <- 20
beta_sparse <- c(3, -2, 1.5, rep(0, pB - 3))     # support vrai = {1,2,3}
true_support <- which(beta_sparse != 0)
sigmaB <- 1
lambda_grid <- 10^seq(-1, 1.6, length.out = 20)
Rb <- 400

metrics <- data.frame(lambda = lambda_grid, exact = 0, tpr = 0, fpr = 0, size = 0)
for (li in seq_along(lambda_grid)) {
  l <- lambda_grid[li]
  ex <- 0; tpr <- 0; fpr <- 0; sz <- 0
  for (r in seq_len(Rb)) {
    Xb <- matrix(rnorm(nB * pB), nB, pB)
    yb <- as.numeric(Xb %*% beta_sparse + rnorm(nB, sd = sigmaB))
    lf <- lasso_fit(Xb, yb, lambda = l, tol = 1e-9)
    sel <- which(lf$beta != 0)
    ex  <- ex + identical(sel, true_support)
    tpr <- tpr + length(intersect(sel, true_support)) / length(true_support)
    fpr <- fpr + length(setdiff(sel, true_support)) / (pB - length(true_support))
    sz  <- sz + length(sel)
  }
  metrics[li, -1] <- c(ex / Rb, tpr / Rb, fpr / Rb, sz / Rb)
}

cat("\n=== (B) Lasso : récupération du support vrai {1,2,3} ===\n")
best <- which.max(metrics$exact)
cat(sprintf("Meilleur lambda (récup. exacte) = %.3f -> %.1f %% des réplications\n",
            metrics$lambda[best], 100 * metrics$exact[best]))
print(round(metrics[c(1, best, nrow(metrics)), ], 3), row.names = FALSE)
cat("Petit lambda : sur-sélection (FPR élevé) ; grand lambda : sous-sélection (TPR chute).\n")

longB <- reshape(metrics[, c("lambda", "tpr", "fpr", "exact")],
                 varying = c("tpr", "fpr", "exact"), v.names = "taux",
                 timevar = "mesure", times = c("vrais positifs (TPR)",
                 "faux positifs (FPR)", "récupération exacte"), direction = "long")
ggB <- ggplot(longB, aes(lambda, taux, colour = mesure)) +
  geom_line(linewidth = 1) + geom_point(size = 1.2) +
  scale_x_log10() +
  labs(title = "Lasso : sélection de variables selon lambda",
       subtitle = "Fenêtre de lambda où le support vrai {1,2,3} est retrouvé exactement",
       x = expression(lambda), y = "taux", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_04_lasso_support.png"), ggB, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_04_ridge_bv.png, mc_04_lasso_support.png\n", sep = "")
