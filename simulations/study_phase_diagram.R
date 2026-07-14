# =============================================================================
# CONTRIBUTION ORIGINALE — Diagramme de phase des estimateurs d'effet causal
# -----------------------------------------------------------------------------
# Question de recherche : dans le modèle partiellement linéaire
#   y = theta·d + X'beta + eps ,  d = X'gamma + v ,
# QUEL estimateur de theta faut-il préférer, SELON LE RÉGIME ? On confronte cinq
# méthodes de la boîte à outils sur une GRILLE de régimes (dimension p/n ×
# densité des contrôles) et l'on produit une CARTE DE DÉCISION : dans chaque
# cellule, la méthode à la plus faible RMSE PARMI CELLES QUI COUVRENT (>= 0.90).
#
# C'est une synthèse originale : non pas « une méthode marche », mais « laquelle,
# quand, et pourquoi ». Chaque chiffre est accompagné de son erreur Monte Carlo.
# =============================================================================

for (f in c("00_linalg", "01_ols", "04_regularisation", "22_debiased_lasso", "mc_tools"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- outils partagés --------------------------------------------------------
lasso_support <- function(Xs, rc) {
  lam <- 1.1 * sd(rc) * sqrt(2 * nrow(Xs) * log(ncol(Xs)))
  which(lasso_fit(Xs, rc, lambda = lam, standardize = FALSE, intercept = FALSE)$beta != 0)
}
resid_on <- function(M, d) { dc <- d - mean(d); Mc <- sweep(M, 2, colMeans(M))
  Mc - outer(dc, as.numeric(crossprod(dc, Mc)) / sum(dc^2)) }
fit_theta <- function(y, d, X, sel) {
  dd <- data.frame(y = y, d = d); if (length(sel)) dd <- cbind(dd, X[, sel, drop = FALSE])
  names(dd) <- c("y", "d", if (length(sel)) paste0("x", sel))
  ols_summary(ols_fit(as.formula(paste("y ~ d",
    if (length(sel)) paste("+", paste(paste0("x", sel), collapse = "+")) else "")), dd))$coefficients["d", ]
}

# --- les cinq estimateurs (renvoient c(estimate, se)) -----------------------
est_ols <- function(y, d, X, Xs) { p <- ncol(X); n <- length(y); if (p >= n - 2) return(c(NA, NA))
  s <- fit_theta(y, d, X, seq_len(p)); c(s$estimate, s$se) }
est_single <- function(y, d, X, Xs) {
  sS <- lasso_support(scale(resid_on(X, d)), lm(y ~ d)$residuals)
  s <- fit_theta(y, d, X, sS); c(s$estimate, s$se) }
est_double <- function(y, d, X, Xs) {
  sU <- union(lasso_support(Xs, y - mean(y)), lasso_support(Xs, d - mean(d)))
  s <- fit_theta(y, d, X, sU); c(s$estimate, s$se) }
est_debias <- function(y, d, X, Xs) { db <- debiased_lasso(cbind(d, X), y, targets = 1); c(db$estimate[1], db$se[1]) }
est_dml <- function(y, d, X, Xs) { n <- length(y); K <- 2L; fold <- sample(rep_len(1:K, n)); yt <- dt <- numeric(n)
  for (k in 1:K) { te <- which(fold == k); tr <- which(fold != k)
    ly <- lasso_support(Xs[tr, , drop = FALSE], y[tr] - mean(y[tr]))
    md <- lasso_support(Xs[tr, , drop = FALSE], d[tr] - mean(d[tr]))
    by <- if (length(ly)) coef(lm(y[tr] ~ X[tr, ly, drop = FALSE])) else mean(y[tr])
    bd <- if (length(md)) coef(lm(d[tr] ~ X[tr, md, drop = FALSE])) else mean(d[tr])
    Xy <- if (length(ly)) cbind(1, X[te, ly, drop = FALSE]) else cbind(rep(1, length(te)))
    Xd <- if (length(md)) cbind(1, X[te, md, drop = FALSE]) else cbind(rep(1, length(te)))
    yt[te] <- y[te] - Xy %*% by; dt[te] <- d[te] - Xd %*% bd }
  th <- sum(dt * yt) / sum(dt^2); psi <- dt * (yt - th * dt)
  c(th, sqrt(mean(psi^2) / (n * mean(dt^2)^2))) }

METHODS <- list(`OLS complet` = est_ols, `sélection simple` = est_single,
                `sélection double` = est_double, `lasso débiaisé` = est_debias, `DML` = est_dml)

gen <- function(n, p, s, theta = 1) { X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(0.4, s), rep(0, p - s)); gamma <- c(rep(0.5, s), rep(0, p - s))
  d <- as.numeric(X %*% gamma) + rnorm(n); y <- theta * d + as.numeric(X %*% beta) + rnorm(n)
  list(X = X, y = y, d = d) }

# --- grille de régimes ------------------------------------------------------
n <- 150; theta0 <- 1; R <- 300; z <- qnorm(0.975)
grid <- expand.grid(p = c(40, 150, 400), s = c(3, 20))     # dimension × densité
grid$regime_p <- factor(c("40" = "p<n", "150" = "p≈n", "400" = "p>n")[as.character(grid$p)],
                        levels = c("p<n", "p≈n", "p>n"))
grid$regime_s <- factor(ifelse(grid$s <= 5, "creux", "dense"), levels = c("creux", "dense"))

rows <- list()
for (g in seq_len(nrow(grid))) {
  p <- grid$p[g]; s <- grid$s[g]
  est <- se <- matrix(NA, R, length(METHODS), dimnames = list(NULL, names(METHODS)))
  for (r in seq_len(R)) {
    d <- gen(n, p, s); Xs <- scale(d$X)
    for (m in seq_along(METHODS)) { v <- METHODS[[m]](d$y, d$d, d$X, Xs); est[r, m] <- v[1]; se[r, m] <- v[2] }
  }
  for (m in names(METHODS)) {
    e <- est[, m]; s_ <- se[, m]; ok <- is.finite(e)
    if (!any(ok)) { bias <- rmse <- cov <- NA } else {
      bias <- mean(e[ok]) - theta0
      rmse <- sqrt(mean((e[ok] - theta0)^2))
      cov <- mean(abs(e[ok] - theta0) <= z * s_[ok])
    }
    rows[[length(rows) + 1L]] <- data.frame(p = p, s = s, regime_p = grid$regime_p[g],
      regime_s = grid$regime_s[g], methode = m, bias = bias, rmse = rmse,
      coverage = cov, mcse_cov = if (is.na(cov)) NA else sqrt(cov * (1 - cov) / sum(ok)))
  }
}
res <- do.call(rbind, rows)

# --- carte de décision : gagnant = RMSE min PARMI couverture >= 0.90 --------
winner <- do.call(rbind, by(res, list(res$p, res$s), function(d) {
  elig <- d[!is.na(d$coverage) & d$coverage >= 0.90, ]
  best <- if (nrow(elig)) elig$methode[which.min(elig$rmse)] else "aucune (couverture<0.90)"
  data.frame(regime_p = d$regime_p[1], regime_s = d$regime_s[1], p = d$p[1], s = d$s[1], gagnant = best)
}))

cat("=== DIAGRAMME DE PHASE : quel estimateur de theta selon le régime ? ===\n")
cat("   (n =", n, ", theta =", theta0, ", R =", R, " ; couverture nominale 0.95)\n\n")
for (g in seq_len(nrow(winner))) {
  sub <- res[res$p == winner$p[g] & res$s == winner$s[g], ]
  cat(sprintf("-- Régime %s / contrôles %s (p=%d, s=%d) --\n",
              winner$regime_p[g], winner$regime_s[g], winner$p[g], winner$s[g]))
  for (i in seq_len(nrow(sub)))
    cat(sprintf("   %-18s biais %+6.3f | RMSE %5.3f | couv %s\n", sub$methode[i],
                ifelse(is.na(sub$bias[i]), NA, sub$bias[i]), sub$rmse[i],
                ifelse(is.na(sub$coverage[i]), "  -  ", sprintf("%.2f", sub$coverage[i]))))
  cat(sprintf("   => GAGNANT : %s\n\n", winner$gagnant[g]))
}

# --- figure : carte de décision ---------------------------------------------
gg <- ggplot(winner, aes(regime_p, regime_s, fill = gagnant)) +
  geom_tile(colour = "white", linewidth = 2) +
  geom_text(aes(label = gagnant), size = 4) +
  labs(title = "Carte de décision : quel estimateur d'effet causal selon le régime ?",
       subtitle = "gagnant = RMSE minimale parmi les méthodes dont l'IC couvre (≥ 0.90)",
       x = "dimension (p relatif à n = 150)", y = "structure des contrôles", fill = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "study_phase_diagram.png"), gg, width = 8.5, height = 5, dpi = 120)

# --- figure : couverture par méthode et régime ------------------------------
res$cellule <- paste0(res$regime_p, "\n", res$regime_s)
g2 <- ggplot(res[!is.na(res$coverage), ], aes(methode, coverage, fill = methode)) +
  geom_col(show.legend = FALSE) + geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_wrap(~ cellule, nrow = 2) + coord_flip(ylim = c(0, 1)) +
  labs(title = "Couverture des IC à 95 % par méthode et régime",
       x = NULL, y = "couverture empirique") + theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "study_phase_coverage.png"), g2, width = 10, height = 6, dpi = 120)

saveRDS(list(res = res, winner = winner), file.path(out_dir, "study_phase_diagram.rds"))
cat("Graphiques -> study_phase_diagram.png, study_phase_coverage.png\n")
