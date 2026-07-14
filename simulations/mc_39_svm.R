# =============================================================================
# Monte Carlo — Module 39 : l'astuce du noyau et l'effet de C
# (1) Sur des donnees non lineairement separables (cercles concentriques), le SVM
#     lineaire echoue la ou le RBF reussit. (2) C regle la parcimonie : petit C =
#     marge large, beaucoup de vecteurs de support ; grand C = ajustement serre.
# =============================================================================

source(file.path("R", "39_svm.R"))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gen_circles <- function(n) {                                # deux cercles concentriques
  r <- c(runif(n/2, 0, 1), runif(n/2, 2, 3)); th <- runif(n, 0, 2*pi)
  X <- cbind(r*cos(th), r*sin(th)); y <- rep(c(1, -1), each = n/2)
  list(X = X, y = y)
}

## (1) lineaire vs RBF --------------------------------------------------------
R <- 30; acc_lin <- acc_rbf <- 0
for (r in seq_len(R)) {
  tr <- gen_circles(200); te <- gen_circles(400)
  acc_lin <- acc_lin + mean(svm_predict(svm_fit(tr$X, tr$y, C = 1, "linear"), te$X) == te$y)
  acc_rbf <- acc_rbf + mean(svm_predict(svm_fit(tr$X, tr$y, C = 1, "rbf", gamma = 1), te$X) == te$y)
}
cat("=== (1) Cercles concentriques (non lineairement separables) ===\n")
cat(sprintf("  SVM lineaire : precision test %.3f  (echoue)\n", acc_lin / R))
cat(sprintf("  SVM RBF      : precision test %.3f  (l'astuce du noyau reussit)\n\n", acc_rbf / R))

## (2) effet de C sur le nombre de vecteurs de support ------------------------
set.seed(1); d <- gen_circles(200); Cs <- c(0.01, 0.1, 1, 10, 100)
cat("=== (2) Effet de C (RBF, gamma=1) : marge vs parcimonie ===\n")
cat(sprintf("%8s %14s\n", "C", "nb vect. support"))
nsv <- sapply(Cs, function(cc) svm_fit(d$X, d$y, C = cc, "rbf", gamma = 1)$n_sv)
for (i in seq_along(Cs)) cat(sprintf("%8.2f %14d\n", Cs[i], nsv[i]))
cat("\n=> Petit C : marge large, beaucoup de vecteurs de support (robuste, lisse).\n")
cat("   Grand C : ajustement serre, moins de support (frontiere plus complexe).\n")

# figure : frontiere RBF sur un tirage
d1 <- gen_circles(300); m <- svm_fit(d1$X, d1$y, C = 1, "rbf", gamma = 1)
gr <- expand.grid(x = seq(-3.2, 3.2, length.out = 120), y = seq(-3.2, 3.2, length.out = 120))
gr$z <- svm_predict(m, as.matrix(gr), decision = TRUE)
gg <- ggplot() +
  geom_raster(data = gr, aes(x, y, fill = z), alpha = 0.6) +
  scale_fill_gradient2(low = "#00798c", mid = "white", high = "#d1495b") +
  geom_contour(data = gr, aes(x, y, z = z), breaks = 0, colour = "black") +
  geom_point(data = data.frame(x = d1$X[,1], y = d1$X[,2], cl = factor(d1$y)),
             aes(x, y, shape = cl), size = 1) +
  labs(title = "SVM à noyau RBF : frontière non linéaire (cercles concentriques)",
       subtitle = "le noyau sépare ce qu'aucun hyperplan ne peut", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "none") + coord_equal()
ggsave(file.path(out_dir, "mc_39_svm.png"), gg, width = 7, height = 6, dpi = 120)
cat("\nGraphique -> mc_39_svm.png\n")
