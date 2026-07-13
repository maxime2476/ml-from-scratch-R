# =============================================================================
# Étude de simulation rigoureuse — convergence et taux
# Complément transversal aux 18 Monte Carlo : on vérifie la CONSISTANCE
# (biais -> 0) et le TAUX sqrt(n) des estimateurs, avec ERREURS MONTE CARLO
# systématiques, et la couverture jugée à son erreur MC près.
# =============================================================================

for (f in c("00_linalg", "01_ols", "03_glm_irls", "mc_tools"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ns <- c(50L, 100L, 200L, 400L, 800L, 1600L); R <- 1000L

# ---- (A) sqrt(n)-consistance de trois estimateurs ---------------------------
sim_ols <- function(n) {                     # coef OLS de x1 (vrai = 2)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(n)
  ols_fit(y ~ x1 + x2, d)$coefficients["x1"]
}
sim_logit <- function(n) {                   # MLE logistique de x1 (vrai = 1)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- rbinom(n, 1, plogis(0.5 + 1 * d$x1 - 0.7 * d$x2))
  glm_irls(y ~ x1 + x2, d, "binomial")$coefficients["x1"]
}
sim_ridge <- function(n) {                   # ridge lambda = sqrt(n) (vrai = 2)
  X <- cbind(rnorm(n), rnorm(n)); y <- X %*% c(2, -1) + rnorm(n)
  lam <- sqrt(n)
  as.numeric(solve(crossprod(X) + lam * diag(2), crossprod(X, y)))[1]
}

studies <- list(
  "OLS (theta=2)"        = list(fn = sim_ols,   truth = 2),
  "Logistique (theta=1)" = list(fn = sim_logit, truth = 1),
  "Ridge lambda=sqrt(n)" = list(fn = sim_ridge, truth = 2))

cat("=== (A) Convergence et taux (erreurs Monte Carlo entre parenthèses) ===\n")
allconv <- list()
for (nm in names(studies)) {
  conv <- convergence_study(studies[[nm]]$fn, ns, R, studies[[nm]]$truth, seed = 1)
  conv$estimateur <- nm; allconv[[nm]] <- conv
  cat("\n--", nm, "-- pente log-log RMSE =", round(rmse_rate(conv), 3), "(cible -0.5)\n")
  tab <- data.frame(
    n = conv$n,
    biais = sprintf("%+.4f (%.4f)", conv$bias, conv$bias_se),
    RMSE = sprintf("%.4f (%.4f)", conv$rmse, conv$rmse_se),
    `sqrt(n)*sd` = round(conv$sqrtn_sd, 3), check.names = FALSE)
  print(tab, row.names = FALSE)
}
convdf <- do.call(rbind, allconv)

ggA <- ggplot(convdf, aes(n, rmse, colour = estimateur)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = rmse - 1.96 * rmse_se, ymax = rmse + 1.96 * rmse_se), width = 0.03) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Convergence : RMSE en fonction de n (échelle log-log)",
       subtitle = "Pente ~ -0.5 = taux sqrt(n) ; barres = erreur Monte Carlo",
       x = "n (log)", y = "RMSE (log)", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_convergence_rmse.png"), ggA, width = 8, height = 5, dpi = 120)

# ---- (B) Normalité asymptotique : sqrt(n)*(theta - theta0) -------------------
cat("\n=== (B) sqrt(n)*(theta_hat - theta0) : loi limite (OLS) ===\n")
nn <- c(50L, 800L); dens <- list()
for (n in nn) {
  s <- sqrt(n) * (vapply(seq_len(2000), function(r) sim_ols(n), numeric(1)) - 2)
  dens[[as.character(n)]] <- s
  cat(sprintf("n=%4d : moyenne=%+.3f  sd=%.3f\n", n, mean(s), sd(s)))
}
dfB <- do.call(rbind, lapply(names(dens), function(nm)
  data.frame(z = dens[[nm]], n = paste0("n = ", nm))))
ggB <- ggplot(dfB, aes(z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "grey80", colour = "white") +
  stat_function(fun = dnorm, args = list(sd = sd(dens[["800"]])), colour = "#d73027", linewidth = 1) +
  facet_wrap(~ n) +
  labs(title = "Normalité asymptotique : sqrt(n)(theta_hat - theta0)",
       subtitle = "La loi limite ne dépend plus de n (même dispersion) — estimateur sqrt(n)-consistant",
       x = expression(sqrt(n)(hat(theta) - theta[0])), y = "densité") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_convergence_clt.png"), ggB, width = 8, height = 4, dpi = 120)

# ---- (C) Couverture jugée à son erreur Monte Carlo --------------------------
cat("\n=== (C) Couverture des IC à 95 % (OLS) avec verdict à l'erreur MC ===\n")
for (n in c(30L, 100L, 500L)) {
  cov <- vapply(seq_len(2000), function(r) {
    d <- data.frame(x1 = rnorm(n)); d$y <- 1 + 2 * d$x1 + rnorm(n)
    fit <- ols_fit(y ~ x1, d); ci <- ols_confint(fit)["x1", ]
    ci[1] <= 2 && 2 <= ci[2]
  }, logical(1))
  cm <- coverage_mc(cov, 0.95)
  cat(sprintf("n=%3d : couverture %.3f ± %.3f  IC=[%.3f,%.3f]  nominal 0.95 %s\n",
              n, cm$coverage, cm$se, cm$ci[1], cm$ci[2],
              if (cm$nominal_ok) "atteint (dans l'erreur MC)" else "NON atteint"))
}

cat("\nGraphiques -> ", out_dir, "/mc_convergence_rmse.png, mc_convergence_clt.png\n", sep = "")
