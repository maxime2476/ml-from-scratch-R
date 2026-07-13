# =============================================================================
# Monte Carlo — Module 17 : couverture des IC bootstrap
# Sur une statistique à loi d'échantillonnage ASYMÉTRIQUE et petit échantillon,
# comparer la couverture des IC : asymptotique (normal), percentile, basique,
# BCa. BCa corrige biais et asymétrie -> meilleure couverture.
#
# DGP : X ~ lognormale ; on estime la moyenne de population theta = exp(mu+sig^2/2).
# =============================================================================

for (f in c("00_linalg", "17_bootstrap", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mu <- 0; sig <- 0.8
theta <- exp(mu + sig^2 / 2)                 # vraie moyenne de la lognormale
types <- c("normal", "percentile", "basic", "bca")
ns <- c(15L, 40L, 100L); R_mc <- 1000L; R_boot <- 1200L

cat("=== Couverture des IC bootstrap à 95 % (moyenne d'une lognormale) ===\n")
cat("theta vrai =", round(theta, 4), " (asymétrie forte en petit n)\n\n")
res <- list()
for (n in ns) {
  covd <- matrix(FALSE, R_mc, length(types), dimnames = list(NULL, types))
  for (r in seq_len(R_mc)) {
    x <- exp(rnorm(n, mu, sig))
    bt <- bootstrap(x, mean, R = R_boot)
    for (ty in types) { ci <- boot_ci(bt, type = ty); covd[r, ty] <- ci[1] <= theta && theta <= ci[2] }
  }
  row <- data.frame(n = n)
  for (ty in types) { cm <- coverage_mc(covd[, ty]); row[[ty]] <- sprintf("%.3f (%.3f)", cm$coverage, cm$se) }
  res[[as.character(n)]] <- list(n = n, cov = colMeans(covd))
  print(row, row.names = FALSE)
}
cat("\n(entre parenthèses : erreur Monte Carlo). BCa se rapproche le plus de 0.95\n",
    "en petit n ; l'IC normal sous-couvre à cause de l'asymétrie.\n", sep = "")

df <- do.call(rbind, lapply(res, function(z)
  data.frame(n = factor(z$n), type = names(z$cov), couv = as.numeric(z$cov))))
df$type <- factor(df$type, levels = types)
gg <- ggplot(df, aes(type, couv, fill = n)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.7, 1)) +
  labs(title = "Bootstrap : couverture des IC à 95 % (statistique asymétrique)",
       subtitle = "BCa restaure la couverture en petit échantillon ; l'IC normal sous-couvre",
       x = NULL, y = "couverture empirique", fill = "n") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_17_bootstrap.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "mc_17_bootstrap.png"), "\n")
