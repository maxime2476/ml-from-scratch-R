# =============================================================================
# Monte Carlo — Module 32 : arbitrage biais-variance de la fenetre + biais de
# bord (Nadaraya-Watson vs local lineaire). La fenetre h gouverne le compromis ;
# la CV le regle. Au BORD, le local lineaire bat nettement le Nadaraya-Watson.
# =============================================================================

for (f in c("32_nonparametric", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

m0 <- function(x) sin(2 * x)
n <- 200; R <- 200

## (1) EQM integree en fonction de la fenetre (arbitrage biais-variance) -----
bws <- seq(0.05, 0.8, by = 0.05); x0 <- seq(-2.5, 2.5, length.out = 60)
mise <- matrix(0, length(bws), 2)                          # colonnes : NW, local lineaire
for (r in seq_len(R)) {
  x <- sort(runif(n, -3, 3)); y <- m0(x) + 0.3 * rnorm(n)
  for (j in seq_along(bws)) {
    mise[j, 1] <- mise[j, 1] + mean((nadaraya_watson(x, y, x0, bws[j]) - m0(x0))^2)
    mise[j, 2] <- mise[j, 2] + mean((local_linear(x, y, x0, bws[j]) - m0(x0))^2, na.rm = TRUE)
  }
}
mise <- mise / R
cat("=== Arbitrage biais-variance de la fenetre (EQM integree) ===\n")
cat(sprintf("  h optimal Nadaraya-Watson : %.2f (EQM %.4f)\n", bws[which.min(mise[,1])], min(mise[,1])))
cat(sprintf("  h optimal local lineaire  : %.2f (EQM %.4f)\n\n", bws[which.min(mise[,2])], min(mise[,2])))

## (2) Biais de bord : DGP LINEAIRE pres du bord (le local lineaire est exact) --
edge <- 2.9; mlin <- function(x) 1.5 * x; bias_nw <- bias_ll <- 0
for (r in seq_len(R)) {
  x <- sort(runif(n, -3, 3)); y <- mlin(x) + 0.3 * rnorm(n)
  bias_nw <- bias_nw + (nadaraya_watson(x, y, edge, 0.5) - mlin(edge))
  bias_ll <- bias_ll + (local_linear(x, y, edge, 0.5) - mlin(edge))
}
cat(sprintf("=== Biais de BORD (DGP lineaire, x0 = %.1f, h = 0.5) ===\n", edge))
cat(sprintf("  Nadaraya-Watson : biais %+.3f  (biais de bord O(h))\n", bias_nw / R))
cat(sprintf("  local lineaire  : biais %+.3f  (quasi nul : reproduit la pente)\n", bias_ll / R))
cat("\n=> La fenetre arbitre biais-variance ; au bord, le local lineaire corrige\n")
cat("   le biais que le Nadaraya-Watson ne peut eviter.\n")

df <- rbind(data.frame(h = bws, mise = mise[, 1], methode = "Nadaraya-Watson"),
            data.frame(h = bws, mise = mise[, 2], methode = "local lineaire"))
gg <- ggplot(df, aes(h, mise, colour = methode)) + geom_line(linewidth = 1) + geom_point() +
  labs(title = "Regression a noyau : EQM integree en fonction de la fenetre",
       subtitle = "courbe en U (biais-variance) ; le local lineaire domine surtout aux bords",
       x = "fenetre h", y = "EQM integree", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_32_bandwidth.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_32_bandwidth.png\n")
