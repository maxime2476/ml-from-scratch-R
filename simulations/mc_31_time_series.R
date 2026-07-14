# =============================================================================
# Monte Carlo — Module 31 : la regression fallacieuse (Granger-Newbold 1974)
# Regresser deux marches aleatoires INDEPENDANTES l'une sur l'autre produit un
# t "significatif" bien plus de 5 % du temps : une relation FALLACIEUSE. Tester
# la stationnarite (ADF) avant de regresser evite le piege. On compare au cas
# de series STATIONNAIRES independantes (ou le t-test se comporte bien).
# =============================================================================

for (f in c("31_time_series", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

T <- 200; R <- 1000; z <- qnorm(0.975)
rej_rw <- rej_st <- 0; adf_rw <- adf_st <- 0
for (r in seq_len(R)) {
  # deux marches aleatoires INDEPENDANTES
  x1 <- cumsum(rnorm(T)); y1 <- cumsum(rnorm(T))
  s <- summary(lm(y1 ~ x1))$coefficients[2, ]
  rej_rw <- rej_rw + (abs(s["t value"]) > z)              # t-test fallacieux
  adf_rw <- adf_rw + (adf_test(y1)$statistic < -3.4)      # ADF rejette la racine unitaire ?

  # deux series STATIONNAIRES independantes (bruit blanc)
  x2 <- rnorm(T); y2 <- rnorm(T)
  s2 <- summary(lm(y2 ~ x2))$coefficients[2, ]
  rej_st <- rej_st + (abs(s2["t value"]) > z)
  adf_st <- adf_st + (adf_test(y2)$statistic < -3.4)
}
cm <- function(k) sprintf("%.3f (%.3f)", k / R, sqrt((k/R)*(1-k/R)/R))
cat(sprintf("=== Regression fallacieuse (T=%d, R=%d, nominal 5 %%) ===\n\n", T, R))
cat(sprintf("  Deux marches aleatoires independantes :\n"))
cat(sprintf("    taux de rejet du t-test (relation fallacieuse) : %s\n", cm(rej_rw)))
cat(sprintf("    ADF conclut a tort a la stationnarite          : %s\n\n", cm(adf_rw)))
cat(sprintf("  Deux series stationnaires independantes :\n"))
cat(sprintf("    taux de rejet du t-test (correct ~5 %%)          : %s\n", cm(rej_st)))
cat(sprintf("    ADF detecte la stationnarite                   : %s\n", cm(adf_st)))
cat("\n=> Sur des series I(1), le t-test rejette massivement (fallacieux). L'ADF,\n")
cat("   applique AVANT, distingue racine unitaire et stationnarite : d'ou la regle\n")
cat("   'tester la stationnarite avant de regresser'.\n")

df <- data.frame(
  serie = rep(c("marches aleatoires I(1)", "stationnaires I(0)"), each = 1),
  rejet = c(rej_rw / R, rej_st / R))
gg <- ggplot(df, aes(serie, rejet, fill = serie)) +
  geom_col(show.legend = FALSE) + geom_hline(yintercept = 0.05, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Regression fallacieuse : le t-test s'effondre sur des series I(1)",
       subtitle = "taux de rejet du t (nominal 5 %, tirete) pour deux series INDEPENDANTES",
       x = NULL, y = "taux de rejet") + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_31_spurious.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_31_spurious.png\n")
