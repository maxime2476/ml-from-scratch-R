# =============================================================================
# Monte Carlo — Module 24 : les trois visages de Var(IC) coïncident
# La fonction d'influence (24.1) unifie l'inférence : la variance ASYMPTOTIQUE
# (sandwich, M2/14), le BOOTSTRAP (M17) et le JACKKNIFE estiment tous la même
# quantité Var(IC)/n. On le vérifie contre la variance d'échantillonnage VRAIE
# (obtenue par simulation), pour la pente d'une régression hétéroscédastique.
# =============================================================================

for (f in c("00_linalg", "01_ols", "02_gls_robust", "17_bootstrap", "24_influence", "mc_tools"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- 200; R <- 1000
slope_hat <- se_if <- se_jack <- se_boot <- numeric(R)
for (r in seq_len(R)) {
  x <- rnorm(n); e <- rnorm(n) * (1 + abs(x))          # hétéroscédasticité
  y <- 1 + 2 * x + e
  X <- cbind(1, x); dd <- data.frame(y = y, x = x)
  slope_hat[r] <- coef(lm(y ~ x))[2]
  se_if[r]   <- sqrt(influence_ols(X, y)$vcov[2, 2])                        # IC = sandwich
  se_jack[r] <- sqrt(jackknife(dd, function(z) coef(lm(y ~ x, z))[2])$var)  # jackknife
  if (r <= 200)                                                            # bootstrap : coûteux
    se_boot[r] <- sd(bootstrap(dd, function(z) coef(lm(y ~ x, z))[2], R = 400)$replicates)
}
truth_sd <- sd(slope_hat)                              # variance d'échantillonnage VRAIE

cat("=== Erreur standard de la pente : Var(IC) sous trois formes ===\n\n")
cat(sprintf("  vérité (sd empirique sur %d simulations) : %.4f\n", R, truth_sd))
cat(sprintf("  fonction d'influence / sandwich (M2)      : %.4f\n", mean(se_if)))
cat(sprintf("  jackknife                                 : %.4f\n", mean(se_jack)))
cat(sprintf("  bootstrap (M17)                           : %.4f\n", mean(se_boot[se_boot > 0])))
cat("\n=> Les trois estimateurs de Var(IC) coïncident avec la vérité : ce sont\n")
cat("   trois façons d'estimer le SECOND MOMENT de la même fonction d'influence.\n")

# --- Figure : normalité asymptotique (24.1) + concordance des SE ------------
zst <- (slope_hat - 2) / mean(se_if)
p1 <- ggplot(data.frame(z = zst), aes(sample = z)) +
  stat_qq(size = 0.6, alpha = 0.5) + stat_qq_line(colour = "#d1495b") +
  labs(title = "Normalité asymptotique via la fonction d'influence (éq. 24.1)",
       subtitle = "QQ-plot de (pente - vraie)/se_IC contre la normale",
       x = "quantiles théoriques", y = "quantiles empiriques") +
  theme_minimal(base_size = 12)

est <- data.frame(methode = c("vérité", "influence/sandwich", "jackknife", "bootstrap"),
                  se = c(truth_sd, mean(se_if), mean(se_jack), mean(se_boot[se_boot > 0])))
est$methode <- factor(est$methode, levels = est$methode)
p2 <- ggplot(est, aes(methode, se, fill = methode)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = truth_sd, linetype = "dashed") +
  labs(title = "Trois estimateurs de Var(IC) = la vérité",
       x = NULL, y = "erreur standard de la pente") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(out_dir, "mc_24_influence.png"), p2, width = 8, height = 5, dpi = 120)
ggsave(file.path(out_dir, "mc_24_normalite.png"), p1, width = 7, height = 5, dpi = 120)
cat("\nGraphiques -> mc_24_influence.png, mc_24_normalite.png\n")
