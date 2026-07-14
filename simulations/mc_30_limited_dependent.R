# =============================================================================
# Monte Carlo — Module 30 : ignorer la structure de la reponse BIAISE l'OLS
# (1) Tobit : l'OLS sur les valeurs censurees rectrecit la pente vers 0 ; le
#     Tobit la recupere. (2) Heckman : sous biais de selection, l'OLS naif est
#     biaise ; la correction de Mills le supprime.
# =============================================================================

for (f in c("03_glm_irls", "30_limited_dependent", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages({ library(ggplot2); has_mass <- requireNamespace("MASS", quietly = TRUE) })
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) Tobit vs OLS sous censure -------------------------------------------
n <- 300; b_true <- 1.0; R <- 400
ols_c <- tob_c <- numeric(R)
for (r in seq_len(R)) {
  x <- rnorm(n); ystar <- 0.5 + b_true * x + rnorm(n); y <- pmax(0, ystar)  # censure a gauche
  ols_c[r] <- coef(lm(y ~ x))[2]
  tob_c[r] <- tobit_fit(y ~ x, data.frame(y = y, x = x), left = 0)$coefficients[2]
}
cat(sprintf("=== (1) Censure : pente vraie = %.1f ===\n", b_true))
cat(sprintf("  OLS (ignore la censure) : %.3f  (biais %+.3f)\n", mean(ols_c), mean(ols_c) - b_true))
cat(sprintf("  Tobit (MLE)             : %.3f  (biais %+.3f)\n\n", mean(tob_c), mean(tob_c) - b_true))

## (2) Heckman vs OLS naif sous biais de selection --------------------------
if (has_mass) {
  ols_h <- heck_h <- numeric(R)
  for (r in seq_len(R)) {
    N <- 800; zs <- rnorm(N); xo <- rnorm(N)
    E <- MASS::mvrnorm(N, c(0, 0), matrix(c(1, 0.7, 0.7, 1), 2))
    d <- as.integer(0.2 + 0.9 * zs + 0.4 * xo + E[, 1] > 0)
    y <- ifelse(d == 1, 1 + 1.0 * xo + E[, 2], NA)                 # pente xo = 1
    dat <- data.frame(d = d, zs = zs, xo = xo, y = y)
    ols_h[r] <- coef(lm(y ~ xo, dat, subset = d == 1))[2]         # OLS naif (selectionnes)
    heck_h[r] <- heckman(d ~ zs + xo, y ~ xo, dat)$beta["xo"]
  }
  cat("=== (2) Biais de selection : pente vraie de xo = 1.0 ===\n")
  cat(sprintf("  OLS naif (selectionnes) : %.3f  (biais %+.3f)\n", mean(ols_h), mean(ols_h) - 1))
  cat(sprintf("  Heckman (correction Mills): %.3f  (biais %+.3f)\n", mean(heck_h), mean(heck_h) - 1))
}
cat("\n=> Modeliser la structure de la reponse (censure, selection) supprime le biais.\n")

df <- data.frame(est = c(ols_c, tob_c), methode = rep(c("OLS (censure ignoree)", "Tobit"), each = R))
gg <- ggplot(df, aes(est, fill = methode)) + geom_density(alpha = 0.5) +
  geom_vline(xintercept = b_true, linetype = "dashed") +
  labs(title = "Censure : l'OLS retrecit la pente, le Tobit la recupere",
       subtitle = paste0("pente vraie = ", b_true, " (tiretee)"), x = "pente estimee", y = "densite", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_30_tobit.png"), gg, width = 8, height = 5, dpi = 120)
cat("Graphique -> mc_30_tobit.png\n")
