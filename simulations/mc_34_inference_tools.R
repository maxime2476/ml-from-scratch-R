# =============================================================================
# Monte Carlo — Module 34 : FWER (Bonferroni) vs FDR (Benjamini-Hochberg)
# On teste m hypotheses dont m0 vraies (H0) et m1 fausses (effet reel). On mesure :
#   FWER : P(au moins un faux positif) ;  FDR : E(faux positifs / total rejets).
# Bonferroni controle le FWER (tres conservateur) ; BH controle le FDR a sa cible
# tout en detectant BEAUCOUP plus de vrais effets (puissance).
# =============================================================================

for (f in c("34_inference_tools", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

m <- 200; m1 <- 40; mu <- 3.2; alpha <- 0.05; R <- 2000
fp_bonf <- fp_bh <- tp_bonf <- tp_bh <- fdr_bonf <- fdr_bh <- 0
for (r in seq_len(R)) {
  z <- c(rnorm(m1, mu), rnorm(m - m1, 0))                   # m1 vrais effets, reste sous H0
  p <- 2 * pnorm(-abs(z)); truth <- c(rep(TRUE, m1), rep(FALSE, m - m1))
  rej_b <- p_adjust_bonferroni(p) < alpha; rej_h <- p_adjust_bh(p) < alpha
  fp_bonf <- fp_bonf + any(rej_b & !truth); fp_bh <- fp_bh + any(rej_h & !truth)   # >=1 faux positif
  tp_bonf <- tp_bonf + sum(rej_b & truth) / m1; tp_bh <- tp_bh + sum(rej_h & truth) / m1
  fdr_bonf <- fdr_bonf + (if (sum(rej_b)) sum(rej_b & !truth) / sum(rej_b) else 0)
  fdr_bh <- fdr_bh + (if (sum(rej_h)) sum(rej_h & !truth) / sum(rej_h) else 0)
}
cat(sprintf("=== Tests multiples (m=%d, dont %d vrais effets, alpha=%.2f, R=%d) ===\n\n", m, m1, alpha, R))
cat(sprintf("%-16s %10s %10s %12s\n", "methode", "FWER", "FDR", "puissance"))
cat(sprintf("%-16s %10.3f %10.3f %12.3f\n", "Bonferroni", fp_bonf/R, fdr_bonf/R, tp_bonf/R))
cat(sprintf("%-16s %10.3f %10.3f %12.3f\n", "Benjamini-Hochberg", fp_bh/R, fdr_bh/R, tp_bh/R))
cat(sprintf("\n=> Bonferroni : FWER <= %.2f (tres conservateur, faible puissance).\n", alpha))
cat(sprintf("   BH : FDR controle a ~%.2f, mais DETECTE bien plus de vrais effets.\n", alpha))

df <- data.frame(
  methode = rep(c("Bonferroni", "Benjamini-Hochberg"), each = 2),
  critere = rep(c("FDR", "puissance"), 2),
  valeur = c(fdr_bonf/R, tp_bonf/R, fdr_bh/R, tp_bh/R))
gg <- ggplot(df, aes(methode, valeur, fill = critere)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  labs(title = "Tests multiples : BH controle le FDR tout en gagnant en puissance",
       subtitle = paste0("FDR cible = ", alpha, " (tirete) ; Bonferroni sacrifie la puissance"),
       x = NULL, y = "taux", fill = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_34_multiple.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_34_multiple.png\n")
