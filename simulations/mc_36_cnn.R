# =============================================================================
# Monte Carlo — Module 36 : un noyau APPREND a detecter un contour vertical par
# descente de gradient (retropropagation de la convolution), et on illustre
# l'ECONOMIE de parametres du partage de poids vs une couche dense.
# =============================================================================

source(file.path("R", "36_cnn.R"))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) Apprendre un detecteur de contour vertical -----------------------------
# Cible : reponse d'un noyau de Sobel vertical. On apprend K a le reproduire.
sobel_v <- array(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), c(3, 3, 1))
gen_img <- function() { X <- matrix(rnorm(12 * 12), 12, 12)
  X[, 1:6] <- X[, 1:6] + 2                                    # bord vertical
  X }
K <- array(rnorm(9) * 0.1, c(3, 3, 1)); lr <- 0.002; loss <- numeric(400)
for (t in seq_len(400)) {
  X <- gen_img(); target <- conv2d(X, sobel_v)$out
  cv <- conv2d(X, K); r <- cv$out - target
  loss[t] <- mean(r^2)
  bw <- conv2d_backward(2 * r / length(r), cv$cache); K <- K - lr * bw$dK
}
cat("=== (1) Apprentissage d'un detecteur de contour ===\n")
cat(sprintf("  perte initiale %.3f -> finale %.4f\n", loss[1], loss[400]))
cat(sprintf("  correlation entre K appris et le noyau de Sobel : %.3f\n\n",
            cor(as.numeric(K), as.numeric(sobel_v))))

## (2) Economie de parametres : conv vs dense ---------------------------------
img <- 28; kk <- 3; Ffilt <- 8; hidden <- 128
p_conv <- kk * kk * Ffilt + Ffilt
p_dense <- img * img * hidden + hidden
cat("=== (2) Parametres pour une image 28x28 ===\n")
cat(sprintf("  couche conv (%d filtres 3x3) : %d parametres\n", Ffilt, p_conv))
cat(sprintf("  couche dense (%d unites)     : %d parametres  (x%.0f)\n\n", hidden, p_dense, p_dense / p_conv))
cat("=> Le partage de poids apprend un detecteur reutilisable partout, avec\n")
cat("   des MILLIERS de fois moins de parametres qu'une couche dense.\n")

gg <- ggplot(data.frame(iter = seq_along(loss), loss = loss), aes(iter, loss)) +
  geom_line(colour = "#00798c", linewidth = 1) + scale_y_log10() +
  labs(title = "CNN : un noyau apprend a detecter un contour vertical",
       subtitle = "perte (echelle log) de la convolution vers la reponse cible (Sobel)",
       x = "iteration", y = "erreur quadratique") + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_36_cnn.png"), gg, width = 8, height = 5, dpi = 120)
cat("Graphique -> mc_36_cnn.png\n")
