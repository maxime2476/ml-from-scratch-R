# =============================================================================
# Monte Carlo — Module 38 : (1) l'attention est une RECHERCHE ASSOCIATIVE dont la
# nettete croit avec l'echelle des similarites ; (2) sans encodage positionnel,
# le Transformer est AVEUGLE a l'ordre (permutation-invariant apres agregation).
# =============================================================================

source(file.path("R", "38_attention.R"))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) Nettete de l'attention vs echelle des cles -----------------------------
# 8 cles orthonormees, une requete alignee sur la cle 3 ; on augmente l'intensite.
d <- 8; K <- diag(d); V <- matrix((1:d) * 5, d, 1)
scales <- c(0.5, 1, 2, 4, 8, 16); wmax <- err <- numeric(length(scales))
for (s in seq_along(scales)) {
  q <- K[3, ] * scales[s]
  a <- attention(matrix(q, 1), K, V)
  wmax[s] <- a$weights[1, 3]; err[s] <- abs(a$out[1, 1] - V[3])   # ecart a la valeur cible V_3=15
}
cat("=== (1) Attention = recherche associative (cle cible = 3, V_3 =", V[3], ") ===\n\n")
cat(sprintf("%8s %14s %14s\n", "echelle", "poids sur cle 3", "erreur"))
for (s in seq_along(scales)) cat(sprintf("%8.1f %14.3f %14.3f\n", scales[s], wmax[s], err[s]))
cat("\n=> Plus la similarite est marquee, plus l'attention se CONCENTRE et recupere\n")
cat("   exactement la bonne valeur (adressage par contenu).\n\n")

## (2) Sans encodage positionnel, l'ordre est invisible -----------------------
# Representation = moyenne des sorties de self-attention (un "pooling" de phrase).
Tn <- 6; dm <- 8
X <- matrix(rnorm(Tn * dm), Tn, dm); perm <- sample(Tn)
repr <- function(Z) colMeans(attention(Z, Z, Z)$out)
diff_noPE <- max(abs(repr(X[perm, ]) - repr(X)))                  # invariant a l'ordre
PE <- positional_encoding(Tn, dm)
diff_PE <- max(abs(repr(X[perm, ] + PE) - repr(X + PE)))          # sensible a l'ordre
cat("=== (2) Sensibilite a l'ORDRE d'une representation de sequence ===\n")
cat(sprintf("  sans encodage positionnel : ecart entre ordres = %.2e  (AVEUGLE)\n", diff_noPE))
cat(sprintf("  avec encodage positionnel : ecart entre ordres = %.3f   (sensible)\n", diff_PE))
cat("=> L'attention seule ignore l'ordre ; l'encodage positionnel le lui donne.\n")

df <- data.frame(echelle = scales, poids = wmax)
gg <- ggplot(df, aes(echelle, poids)) + geom_line(colour = "#00798c", linewidth = 1) + geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") + scale_x_log10() +
  labs(title = "Attention : la nettete croit avec l'intensite de la similarite",
       subtitle = "poids d'attention sur la cle cible (1 = recuperation parfaite)",
       x = "echelle de la requete (log)", y = "poids sur la bonne cle") + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_38_attention.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_38_attention.png\n")
