# =============================================================================
# ÉTUDE MÉTHODOLOGIQUE — Pourquoi la « double sélection » ? (Belloni-Chernozhukov-
# Hansen 2014). Choisir les variables de contrôle par le lasso est tentant, mais
# la façon de le faire change TOUT pour l'inférence causale.
# -----------------------------------------------------------------------------
# Modèle : Y = theta·D + X'beta + eps ,  D = X'gamma + v  (p grand, effets creux).
# On veut theta. Un confondeur X1 prédit FORTEMENT le traitement (gamma_1 grand)
# mais FAIBLEMENT le résultat une fois D contrôlé (beta_1 petit).
#
#  * SÉLECTION SIMPLE : on sélectionne les contrôles par leur association au
#    résultat (en présence de D), puis OLS. -> X1 a un signal résiduel faible,
#    le lasso le RATE ; l'omettre biaise theta (variable omise) et effondre la
#    couverture.
#  * SÉLECTION DOUBLE (BCH) : on prend l'UNION des variables sélectionnées dans
#    l'équation du résultat ET dans celle du traitement, puis OLS. -> X1 est
#    rattrapé par l'équation du traitement (gamma_1 grand) ; theta est valide.
#
# HYPOTHÈSE. Plus le confondeur pèse sur le traitement (gamma_1 croît), plus la
# sélection simple est biaisée et sa couverture s'effondre ; la double sélection
# reste centrée et couvre ~0.95.
# =============================================================================

for (f in c("00_linalg", "01_ols", "04_regularisation", "mc_tools"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- 250; p <- 100; theta0 <- 1.0; beta1 <- 0.15

# support lasso avec pénalité théorique calée sur le bruit du résidu
lasso_support <- function(Xs, rc) {
  lam <- 1.1 * sd(rc) * sqrt(2 * nrow(Xs) * log(ncol(Xs)))
  which(lasso_fit(Xs, rc, lambda = lam, standardize = FALSE, intercept = FALSE)$beta != 0)
}
# résidualisation vectorisée de chaque colonne de M sur d (partialling-out, FWL)
resid_on <- function(M, d) {
  dc <- d - mean(d); Mc <- sweep(M, 2, colMeans(M)); Mc - outer(dc, as.numeric(crossprod(dc, Mc)) / sum(dc^2))
}
fit_theta <- function(y, d, X, sel) {
  dd <- data.frame(y = y, d = d); if (length(sel)) dd <- cbind(dd, X[, sel, drop = FALSE])
  names(dd) <- c("y", "d", if (length(sel)) paste0("x", sel))
  fm <- as.formula(paste("y ~ d", if (length(sel)) paste("+", paste(paste0("x", sel), collapse = "+")) else ""))
  ols_summary(ols_fit(fm, dd))$coefficients["d", ]
}

gammas <- c(0.5, 1.5, 2.5); R <- 300; z <- qnorm(0.975)
rows <- list()
for (g1 in gammas) {
  bss <- bds <- numeric(R); css <- cds <- 0
  for (r in seq_len(R)) {
    X <- matrix(rnorm(n * p), n, p)
    beta <- c(beta1, rep(0, p - 1)); gamma <- c(g1, rep(0, p - 1))
    d <- as.numeric(X %*% gamma) + rnorm(n); y <- theta0 * d + as.numeric(X %*% beta) + rnorm(n)
    Xs <- scale(X)
    # simple : sélectionner en contrôlant pour D (résidus sur D)
    sS <- lasso_support(scale(resid_on(X, d)), as.numeric(lm(y ~ d)$residuals))
    # double : union des supports Y~X et D~X
    sU <- union(lasso_support(Xs, y - mean(y)), lasso_support(Xs, d - mean(d)))
    es <- fit_theta(y, d, X, sS); ed <- fit_theta(y, d, X, sU)
    bss[r] <- es$estimate; css <- css + (abs(es$estimate - theta0) <= z * es$se)
    bds[r] <- ed$estimate; cds <- cds + (abs(ed$estimate - theta0) <= z * ed$se)
  }
  rows[[length(rows) + 1L]] <- data.frame(
    gamma1 = g1, methode = c("sélection simple", "sélection double (BCH)"),
    biais = c(mean(bss) - theta0, mean(bds) - theta0),
    couverture = c(css / R, cds / R))
}
tab <- do.call(rbind, rows)

cat(sprintf("=== Double sélection (PLM haute dim : n=%d, p=%d, theta0=%.1f, R=%d) ===\n\n",
            n, p, theta0, R))
cat(sprintf("%8s %-24s %9s %11s\n", "gamma1", "méthode", "biais", "couv. 95%"))
for (i in seq_len(nrow(tab)))
  cat(sprintf("%8.1f %-24s %9.3f %11.3f\n", tab$gamma1[i], tab$methode[i], tab$biais[i], tab$couverture[i]))
cat("\n=> Quand le confondeur pèse sur le traitement (gamma1 croît), la sélection\n")
cat("   SIMPLE le rate : biais croissant, couverture qui s'effondre. La sélection\n")
cat("   DOUBLE le rattrape par l'équation du traitement : biais ~0, couverture ~0.95.\n")
cat("   Leçon : sélectionner les contrôles depuis les DEUX équations.\n")

gg <- ggplot(tab, aes(factor(gamma1), couverture, fill = methode)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Double sélection : pourquoi choisir les contrôles depuis les deux équations",
       subtitle = "Couverture des IC à 95 % ; la sélection simple s'effondre quand le confondeur pèse sur D",
       x = expression(gamma[1]~"(force du confondeur sur le traitement)"),
       y = "couverture empirique", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "study_double_selection.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "study_double_selection.png"), "\n")
saveRDS(tab, file.path(out_dir, "study_double_selection.rds"))
