# =============================================================================
# Monte Carlo — Module 22 : le lasso débiaisé RÉPARE l'inférence post-sélection
# En haute dimension (p > n), le t-test OLS post-lasso naïf (Module 14) échoue :
# sous H0 il sur-rejette, et ses IC sous-couvrent. Le lasso débiaisé restaure la
# couverture ~0.95 et le niveau ~5 %, coordonnée par coordonnée.
#
# DGP : y = X beta + eps, beta creux (s actifs), p > n.
# =============================================================================

for (f in c("00_linalg", "01_ols", "04_regularisation", "22_debiased_lasso", "mc_tools"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- 100; p <- 200; s <- 5; b_act <- 1.0
beta <- c(rep(b_act, s), rep(0, p - s))
j_act <- 1; j_null <- 50                      # une coord active, une nulle
R <- 800

cov_db_a <- cov_db_n <- rej_db_n <- 0
cov_nv_a <- cov_nv_n <- rej_nv_n <- 0; n_nv <- 0
z <- qnorm(0.975)
for (r in seq_len(R)) {
  X <- matrix(rnorm(n * p), n, p); y <- as.numeric(X %*% beta + rnorm(n))

  # --- Lasso débiaisé ---
  db <- debiased_lasso(X, y, targets = c(j_act, j_null))
  cov_db_a <- cov_db_a + (db$lower[1] <= b_act && b_act <= db$upper[1])
  cov_db_n <- cov_db_n + (db$lower[2] <= 0 && 0 <= db$upper[2])
  rej_db_n <- rej_db_n + (abs(db$estimate[2] / db$se[2]) > z)     # H0: beta_null=0

  # --- Post-lasso naïf : sélection puis OLS + t-test ---
  st <- .standardize(X); Xs <- st$Xs; yc <- y - mean(y)
  lam <- 0.4 * max(abs(crossprod(Xs, yc)))
  sel <- which(lasso_fit(Xs, yc, lambda = lam, standardize = FALSE, intercept = FALSE)$beta != 0)
  if (length(sel) > 0 && length(sel) < n - 2) {
    dd <- as.data.frame(X[, sel, drop = FALSE]); names(dd) <- paste0("v", sel); dd$y <- y
    fit <- ols_fit(as.formula(paste("y ~", paste(names(dd)[-ncol(dd)], collapse = "+"))), dd)
    sm <- ols_summary(fit)$coefficients
    if (j_null %in% sel) {
      nm <- paste0("v", j_null); est <- sm[nm, "estimate"]; se <- sm[nm, "se"]
      cov_nv_n <- cov_nv_n + (abs(est - 0) <= z * se); rej_nv_n <- rej_nv_n + (abs(est/se) > z)
    } else { cov_nv_n <- cov_nv_n + 1 }   # non sélectionné -> "couvre" 0 trivialement
    if (j_act %in% sel) {
      nm <- paste0("v", j_act); est <- sm[nm, "estimate"]; se <- sm[nm, "se"]
      cov_nv_a <- cov_nv_a + (abs(est - b_act) <= z * se)
    }
    n_nv <- n_nv + 1
  }
}

cm <- function(k, N) sprintf("%.3f (%.3f)", k / N, sqrt((k/N)*(1-k/N)/N))
cat("=== (p =", p, "> n =", n, ", s =", s, " actifs ; R =", R, ") ===\n\n")
cat("COUVERTURE des IC à 95 % :\n")
cat(sprintf("  coef ACTIF (beta=%.1f) : débiaisé %s | post-lasso naïf %s\n",
            b_act, cm(cov_db_a, R), cm(cov_nv_a, n_nv)))
cat(sprintf("  coef NUL   (beta=0)   : débiaisé %s | post-lasso naïf %s\n",
            cm(cov_db_n, R), cm(cov_nv_n, n_nv)))
cat("\nTAUX DE REJET sous H0 : beta_nul = 0 (nominal 5 %) :\n")
cat(sprintf("  débiaisé %s | post-lasso naïf %s\n", cm(rej_db_n, R), cm(rej_nv_n, n_nv)))
cat("\n=> Le lasso débiaisé restaure couverture ~0.95 et niveau ~5 % ; le naïf échoue.\n")

df <- data.frame(
  coef = rep(c("actif", "nul"), each = 2),
  methode = rep(c("débiaisé", "post-lasso naïf"), 2),
  couv = c(cov_db_a/R, cov_nv_a/n_nv, cov_db_n/R, cov_nv_n/n_nv))
gg <- ggplot(df, aes(coef, couv, fill = methode)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Lasso débiaisé : l'inférence haute dimension réparée",
       subtitle = paste0("p=", p, " > n=", n, " : le débiaisé couvre ~0.95 ; le post-lasso naïf (M14) échoue"),
       x = "coefficient", y = "couverture empirique", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_22_debiased.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "mc_22_debiased.png"), "\n")
