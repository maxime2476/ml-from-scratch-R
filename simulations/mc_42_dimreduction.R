# =============================================================================
# Monte Carlo — Module 42 : (1) la PCA a noyau LINEARISE des cercles concentriques
# que la PCA lineaire ne peut separer ; (2) l'ICA recupere des sources
# independantes la ou la PCA ne fait que decorreler.
# =============================================================================

for (f in c("42_dimreduction", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) kernel PCA vs PCA lineaire sur cercles concentriques --------------------
R <- 40; sep_lin <- sep_ker <- 0
for (r in seq_len(R)) {
  n <- 200; th <- runif(n, 0, 2*pi); rad <- c(rep(1, n/2), rep(4, n/2)) + rnorm(n, 0, 0.15)
  X <- cbind(rad * cos(th), rad * sin(th)); y <- rep(1:2, each = n/2)
  sep_lin <- sep_lin + abs(cor(prcomp(X)$x[, 1], y))          # PCA lineaire : 1re composante
  sep_ker <- sep_ker + abs(cor(kernel_pca(X, 1, gamma = 0.5)$proj[, 1], y))  # kernel PCA
}
cat("=== (1) Separer deux cercles concentriques (|correlation| composante 1) ===\n\n")
cat(sprintf("  PCA lineaire  : %.3f  (aucune direction lineaire ne separe les rayons)\n", sep_lin / R))
cat(sprintf("  PCA a noyau   : %.3f  (deplie les cercles en une coordonnee radiale)\n\n", sep_ker / R))

## (2) ICA vs PCA : recuperer des sources independantes ------------------------
R <- 30; cor_ica <- cor_pca <- 0
for (r in seq_len(R)) {
  n <- 800; s1 <- sin((1:n)/19); s2 <- 2*(((1:n)*7 %% 100)/100 - 0.5)
  S <- scale(cbind(s1, s2))                                    # variances EGALES -> PCA ambigue
  ang <- pi/4; A <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), 2)  # rotation 45
  X <- S %*% t(A)
  match_src <- function(Z) mean(sapply(1:2, function(k) max(abs(cor(Z, S[, k])))))
  cor_ica <- cor_ica + match_src(ica_fastica(X, 2)$S)
  cor_pca <- cor_pca + match_src(prcomp(X)$x)
}
cat("=== (2) Recuperer des sources INDEPENDANTES (|correlation| aux vraies sources) ===\n")
cat(sprintf("  PCA (decorrele seulement) : %.3f\n", cor_pca / R))
cat(sprintf("  ICA (independance)        : %.3f\n", cor_ica / R))
cat("\n=> La PCA a noyau capte la structure courbe ; l'ICA separe les sources la ou\n")
cat("   la PCA, limitee aux moments d'ordre 2, echoue.\n")

# figure : kernel PCA sur cercles
n <- 300; th <- runif(n, 0, 2*pi); rad <- c(rep(1, n/2), rep(4, n/2)) + rnorm(n, 0, 0.12)
X <- cbind(rad*cos(th), rad*sin(th)); y <- factor(rep(1:2, each = n/2))
kp <- kernel_pca(X, 2, 0.5)
df <- rbind(data.frame(x = X[,1], y2 = X[,2], cl = y, vue = "espace original"),
            data.frame(x = kp$proj[,1], y2 = kp$proj[,2], cl = y, vue = "PCA a noyau"))
gg <- ggplot(df, aes(x, y2, colour = cl)) + geom_point(size = 1) + facet_wrap(~vue, scales = "free") +
  labs(title = "PCA a noyau : deplie des cercles concentriques en clusters lineairement separables",
       x = NULL, y = NULL, colour = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_42_kpca.png"), gg, width = 9, height = 4.5, dpi = 120)
cat("\nGraphique -> mc_42_kpca.png\n")
