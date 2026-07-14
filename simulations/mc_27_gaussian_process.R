# =============================================================================
# Monte Carlo — Module 27 : le GP quantifie l'incertitude et choisit l'échelle
# (1) Couverture des intervalles a posteriori à 95 % quand le GP est le vrai DGP.
# (2) La vraisemblance marginale (rasoir d'Occam) sélectionne la bonne échelle l
#     — alternative principielle à la validation croisée (M6).
# =============================================================================

for (f in c("00_linalg", "27_gaussian_process")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) Couverture : tirer f d'un GP, observer avec bruit, vérifier les IC ----
l0 <- 1; sf0 <- 1; sn0 <- 0.2; n <- 40
sample_gp <- function(Xall) {                              # tire une fonction du GP a priori
  K <- rbf_kernel(Xall, Xall, l0, sf0^2) + 1e-8 * diag(nrow(Xall))
  as.numeric(t(chol(K)) %*% rnorm(nrow(Xall)))
}
R <- 500; z <- qnorm(0.975); cov_f <- 0; ntest <- 0
for (r in seq_len(R)) {
  Xall <- matrix(sort(runif(n + 20, 0, 10)), n + 20, 1)
  f <- sample_gp(Xall); tr <- 1:n; te <- (n + 1):(n + 20)
  y <- f[tr] + sn0 * rnorm(n)
  pr <- gp_predict(gp_fit(Xall[tr, , drop = FALSE], y, l0, sf0, sn0), Xall[te, , drop = FALSE])
  cov_f <- cov_f + sum(abs(f[te] - pr$mean) <= z * pr$sd); ntest <- ntest + length(te)
}
cat("=== (1) Couverture des IC a posteriori du GP (95 % nominal) ===\n")
cat(sprintf("  couverture de f(x*) : %.3f  (sur %d prédictions)\n\n", cov_f / ntest, ntest))

## (2) Sélection de l'échelle par vraisemblance marginale --------------------
cat("=== (2) La vraisemblance marginale sélectionne l'échelle (vraie l0 = 1) ===\n\n")
set.seed(7); Xtr <- matrix(sort(runif(50, 0, 10)), 50, 1)
ytr <- sample_gp(Xtr) + sn0 * rnorm(50)
ls <- seq(0.2, 4, by = 0.1)
ll <- sapply(ls, function(l) gp_fit(Xtr, ytr, l, sf0, sn0)$loglik)
l_hat <- ls[which.max(ll)]
opt <- gp_optimize(Xtr, ytr)
cat(sprintf("  échelle sur grille (argmax loglik) : %.2f\n", l_hat))
cat(sprintf("  échelle par optim (L-BFGS-B)       : %.2f  (vraie ~ %.1f)\n", opt$lengthscale, l0))
cat("  => l'automatisme bayésien retrouve l'échelle, sans validation croisée.\n")

gg <- ggplot(data.frame(l = ls, ll = ll), aes(l, ll)) +
  geom_line(linewidth = 1, colour = "#00798c") +
  geom_vline(xintercept = l0, linetype = "dashed") +
  geom_point(data = data.frame(l = l_hat, ll = max(ll)), colour = "#d1495b", size = 3) +
  labs(title = "Vraisemblance marginale : le rasoir d'Occam sélectionne l'échelle",
       subtitle = paste0("maximum à l = ", l_hat, " (vraie échelle = ", l0, ", tiretée)"),
       x = expression("échelle "*ell), y = "log-vraisemblance marginale") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_27_marginal_likelihood.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_27_marginal_likelihood.png\n")
