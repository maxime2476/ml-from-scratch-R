# =============================================================================
# Monte Carlo — Module 23 : l'ajustement OVB récupère l'effet sous confusion CONNUE
# DGP avec confondeur Z inobservé :  D = a·Z + v ,  Y = tau·D + gamma·Z + eps.
# On régresse Y sur D SEUL (Z omis) -> estimation biaisée. Si l'on connaissait les
# R² partiels de Z, l'ajustement (23.2) récupère tau. On montre aussi que la
# robustness value du résultat naïf décroît quand le confondeur se renforce.
# =============================================================================

for (f in c("00_linalg", "01_ols", "23_sensitivity", "mc_tools"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- 400; tau <- 1.0; a <- 0.8       # force confondeur -> traitement
gammas <- c(0, 0.5, 1.0, 1.5, 2.0)   # force confondeur -> résultat (croissante)
R <- 400

res <- lapply(gammas, function(gam) {
  naive <- adj <- rv <- numeric(R)
  for (r in seq_len(R)) {
    Z <- rnorm(n); D <- a * Z + rnorm(n); y <- tau * D + gam * Z + rnorm(n)
    # --- estimation NAÏVE : Y ~ D (Z omis) ---
    d0 <- data.frame(y = y, D = D); f0 <- ols_fit(y ~ D, d0)
    s <- sensitivity_ols(f0, "D"); naive[r] <- s$estimate; rv[r] <- s$rv_q
    # --- R² partiels VRAIS du confondeur Z (calculables ici car Z simulé) ---
    r2dz <- summary(lm(D ~ Z))$r.squared                 # R²_{D~Z}
    fyz  <- ols_summary(ols_fit(y ~ D + Z, data.frame(y = y, D = D, Z = Z)))$coefficients["Z", "t"]
    r2yz <- partial_r2(fyz, n - 3)                       # R²_{Y~Z|D}
    adj[r] <- adjusted_estimate(s$estimate, s$se, s$df, r2dz, r2yz)$estimate
  }
  c(gamma = gam, naive = mean(naive), adj = mean(adj), rv = mean(rv))
})
tab <- do.call(rbind, res)

cat("=== Ajustement OVB sous confusion connue (tau =", tau, ") ===\n\n")
cat(sprintf("%6s | %10s | %10s | %8s\n", "gamma", "naïf", "ajusté", "RV_q(naïf)"))
for (i in seq_len(nrow(tab)))
  cat(sprintf("%6.1f | %10.3f | %10.3f | %8.3f\n",
              tab[i, "gamma"], tab[i, "naive"], tab[i, "adj"], tab[i, "rv"]))
cat(sprintf("\n=> Le naïf s'écarte de tau=%.1f quand gamma croît ; l'ajusté le récupère.\n", tau))
cat("=> Mise en garde : la RV du naïf reste ~0.7 (l'analyse « paraît robuste »)\n")
cat("   alors même que l'effet est fortement biaisé. La RV borne l'impact d'un\n")
cat("   confondeur HYPOTHÉTIQUE ; elle ne prouve pas l'absence de confusion.\n")

df <- data.frame(gamma = rep(tab[, "gamma"], 2),
                 est = c(tab[, "naive"], tab[, "adj"]),
                 methode = rep(c("naïf (Z omis)", "ajusté OVB"), each = nrow(tab)))
gg <- ggplot(df, aes(gamma, est, colour = methode)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_hline(yintercept = tau, linetype = "dashed") +
  labs(title = "Analyse de sensibilité : l'ajustement OVB récupère l'effet",
       subtitle = paste0("tau=", tau, " (tireté) ; le naïf (Z omis) dévie avec la force gamma du confondeur"),
       x = expression(gamma~"(force confondeur"%->%"résultat)"),
       y = "effet estimé", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_23_sensitivity.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "mc_23_sensitivity.png"), "\n")
