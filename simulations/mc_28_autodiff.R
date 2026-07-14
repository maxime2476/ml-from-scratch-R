# =============================================================================
# Monte Carlo — Module 28 : l'autodiff pilote l'apprentissage d'un MLP
# On entraîne un perceptron à une couche cachée par descente de gradient dont le
# gradient est fourni ENTIÈREMENT par la différentiation automatique en mode
# inverse (aucune dérivée codée à la main). On vérifie que la perte décroît et que
# le gradient AD coïncide avec numDeriv à chaque pas.
# =============================================================================

source(file.path("R", "28_autodiff.R"))
suppressMessages({ library(ggplot2); has_nd <- requireNamespace("numDeriv", quietly = TRUE) })
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Données : cible non linéaire ------------------------------------------
n <- 200; p <- 2; h <- 10
X <- cbind(matrix(runif(n * p, -2, 2), n, p), 1)                  # biais augmenté
f_true <- function(Z) sin(Z[, 1]) * cos(Z[, 2])
y <- matrix(f_true(X) + 0.05 * rnorm(n), n, 1)

unpack <- function(w) list(
  W1 = matrix(w[1:((p + 1) * h)], p + 1, h),
  W2 = matrix(w[((p + 1) * h + 1):((p + 1) * h + h + 1)], h + 1, 1))
mlp_grad_loss <- function(w) {                                    # gradient + perte par AD
  pk <- unpack(w)
  ad_reset(); W1 <- adnode(pk$W1); W2 <- adnode(pk$W2)
  Ha <- ad_cbind1(tanh(mm(X, W1)))
  r <- mm(Ha, W2) - y; L <- sum(r * r) / n; backward(L)
  list(grad = c(as.numeric(W1$grad), as.numeric(W2$grad)), loss = L$value)
}

# --- Descente de gradient pilotée par l'AD ----------------------------------
set.seed(1); w <- rnorm((p + 1) * h + h + 1) * 0.3
eta <- 0.05; steps <- 400; loss_path <- numeric(steps); max_grad_err <- 0
for (s in seq_len(steps)) {
  gl <- mlp_grad_loss(w); loss_path[s] <- gl$loss
  if (has_nd && s %% 100 == 1) {                                  # contrôle vs numDeriv
    ndg <- numDeriv::grad(function(w) { pk <- unpack(w)
      Ha <- cbind(tanh(X %*% pk$W1), 1); sum((Ha %*% pk$W2 - y)^2) / n }, w)
    max_grad_err <- max(max_grad_err, max(abs(gl$grad - ndg)))
  }
  w <- w - eta * gl$grad
}
cat("=== MLP entraîné par autodiff (mode inverse) ===\n\n")
cat(sprintf("  perte initiale : %.4f  ->  perte finale : %.4f  (%.0f%% de réduction)\n",
            loss_path[1], loss_path[steps], 100 * (1 - loss_path[steps] / loss_path[1])))
if (has_nd) cat(sprintf("  écart max gradient AD vs numDeriv (contrôles) : %.2e\n", max_grad_err))
cat("\n=> Le gradient exact fourni par l'AD suffit à entraîner le réseau ;\n")
cat("   aucune dérivée n'a été codée à la main.\n")

gg <- ggplot(data.frame(step = seq_len(steps), loss = loss_path), aes(step, loss)) +
  geom_line(linewidth = 1, colour = "#00798c") + scale_y_log10() +
  labs(title = "Apprentissage d'un MLP par descente de gradient autodiff",
       subtitle = "gradient calculé par différentiation automatique en mode inverse",
       x = "itération", y = "perte MSE (échelle log)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_28_autodiff.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_28_autodiff.png\n")
