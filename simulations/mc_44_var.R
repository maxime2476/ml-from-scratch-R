# =============================================================================
# Monte Carlo — Module 44 : le test de Granger detecte la VRAIE direction de
# causalite (taille sous H0, puissance sous H1), et les IRF montrent la
# propagation puis la decroissance d'un choc.
# =============================================================================

for (f in c("44_var", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) Granger : taille et puissance -----------------------------------------
gen <- function(n, beta) { Y <- matrix(0, n, 2)          # y2 -> y1 avec intensite beta
  for (t in 2:n) { Y[t, 1] <- 0.3 * Y[t-1, 1] + beta * Y[t-1, 2] + rnorm(1)
                   Y[t, 2] <- 0.4 * Y[t-1, 2] + rnorm(1) }
  Y }
R <- 1000; n <- 200; a <- 0.05
rej_H0 <- rej_H1 <- 0                                     # H0 : beta=0 (y2 ne cause pas y1)
for (r in seq_len(R)) {
  m0 <- var_fit(gen(n, 0), 1); m1 <- var_fit(gen(n, 0.5), 1)
  rej_H0 <- rej_H0 + (granger_test(m0, cause = 2, effect = 1)$p_value < a)
  rej_H1 <- rej_H1 + (granger_test(m1, cause = 2, effect = 1)$p_value < a)
}
cat(sprintf("=== (1) Test de Granger (n=%d, R=%d, alpha=%.2f) ===\n\n", n, R, a))
cat(sprintf("  taille (H0 : y2 ne cause pas y1)  : %.3f  (~ nominal)\n", rej_H0 / R))
cat(sprintf("  puissance (H1 : y2 cause y1)      : %.3f\n\n", rej_H1 / R))

## (2) IRF : propagation d'un choc -------------------------------------------
set.seed(1); Yd <- gen(400, 0.6); m <- var_fit(Yd, 1); ih <- var_irf(m, 15)
cat("=== (2) Reponse a un choc unitaire en y2 ===\n")
cat("  horizon :", 0:15, "\n")
cat("  y1 :", round(ih[, 1, 2], 2), "\n")
cat("  y2 :", round(ih[, 2, 2], 2), "\n")
cat("=> Le choc en y2 se propage a y1 (causalite), puis les deux reponses\n")
cat("   decroissent vers 0 (stationnarite).\n")

df <- data.frame(h = rep(0:15, 2), irf = c(ih[, 1, 2], ih[, 2, 2]),
                 reponse = rep(c("y1", "y2"), each = 16))
gg <- ggplot(df, aes(h, irf, colour = reponse)) + geom_line(linewidth = 1) + geom_point() +
  geom_hline(yintercept = 0, colour = "grey60") +
  labs(title = "Reponse impulsionnelle : propagation d'un choc en y2",
       subtitle = "le choc atteint y1 (causalite de Granger) puis s'amortit",
       x = "horizon (periodes)", y = "reponse", colour = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_44_irf.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_44_irf.png\n")
