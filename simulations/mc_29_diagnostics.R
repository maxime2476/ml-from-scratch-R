# =============================================================================
# Monte Carlo — Module 29 : taille et puissance des tests de diagnostic
# Sous H0 (modele correct), chaque test doit rejeter au taux nominal (~5 %) :
# c'est la TAILLE. Sous l'alternative (heteroscedasticite / autocorrelation /
# endogeneite), il doit rejeter souvent : c'est la PUISSANCE.
# =============================================================================

for (f in c("29_diagnostics", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- 150; R <- 1000; a <- 0.05
rej <- function(pv) mean(pv < a)

# --- Breusch-Pagan : taille (homosced.) vs puissance (heterosced.) ---------
bp_h0 <- bp_h1 <- numeric(R)
for (r in seq_len(R)) {
  x1 <- rnorm(n); x2 <- rnorm(n)
  y0 <- 1 + x1 - x2 + rnorm(n)                     # homoscedastique
  y1 <- 1 + x1 - x2 + rnorm(n) * exp(0.7 * x1)     # heteroscedastique (forte)
  bp_h0[r] <- bp_test(y0 ~ x1 + x2, data.frame(y0, x1, x2))$p_value
  bp_h1[r] <- bp_test(y1 ~ x1 + x2, data.frame(y1, x1, x2))$p_value
}

# --- Breusch-Godfrey : taille (i.i.d.) vs puissance (AR(1)) -----------------
bg_h0 <- bg_h1 <- numeric(R)
for (r in seq_len(R)) {
  x1 <- rnorm(n); e0 <- rnorm(n); e1 <- as.numeric(filter(rnorm(n), 0.6, "recursive"))
  y0 <- 1 + x1 + e0; y1 <- 1 + x1 + e1
  bg_h0[r] <- bg_test(y0 ~ x1, data.frame(y0, x1), order = 1)$p_value
  bg_h1[r] <- bg_test(y1 ~ x1, data.frame(y1, x1), order = 1)$p_value
}

# --- Durbin-Wu-Hausman : taille (exogene) vs puissance (endogene) -----------
dwh_h0 <- dwh_h1 <- numeric(R)
for (r in seq_len(R)) {
  z <- rnorm(n); x1 <- rnorm(n); v <- rnorm(n); u <- rnorm(n)
  xexo <- 0.7 * z + 0.4 * x1 + v                    # exogene : v ind. de u
  xend <- 0.7 * z + 0.4 * x1 + v; yend <- 1 + 2 * xend + x1 + (u + 0.8 * v)  # endogene: corr(v,err)
  yexo <- 1 + 2 * xexo + x1 + u
  Xo <- cbind(1, xexo, x1); Xe <- cbind(1, xend, x1); Z <- cbind(1, z, x1)
  dwh_h0[r] <- dwh_test(yexo, Xo, Z, endog = 2)$p_value
  dwh_h1[r] <- dwh_test(yend, Xe, Z, endog = 2)$p_value
}

tab <- data.frame(
  test = c("Breusch-Pagan", "Breusch-Godfrey", "Durbin-Wu-Hausman"),
  taille = c(rej(bp_h0), rej(bg_h0), rej(dwh_h0)),
  puissance = c(rej(bp_h1), rej(bg_h1), rej(dwh_h1)))

cat(sprintf("=== Taille et puissance des tests (n=%d, R=%d, alpha=%.2f) ===\n\n", n, R, a))
cat(sprintf("%-20s %8s %10s\n", "test", "taille", "puissance"))
for (i in seq_len(nrow(tab)))
  cat(sprintf("%-20s %8.3f %10.3f\n", tab$test[i], tab$taille[i], tab$puissance[i]))
cat(sprintf("\n(erreur MC binomiale ~ %.3f) => taille ~ 0.05, puissance elevee sous l'alternative.\n",
            sqrt(0.05 * 0.95 / R)))

df <- rbind(data.frame(test = tab$test, taux = tab$taille, regime = "taille (H0)"),
            data.frame(test = tab$test, taux = tab$puissance, regime = "puissance (H1)"))
gg <- ggplot(df, aes(test, taux, fill = regime)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = a, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Tests de diagnostic : taille nominale sous H0, puissance sous H1",
       subtitle = "la taille colle à 0.05 (tireté) ; la puissance détecte le défaut",
       x = NULL, y = "taux de rejet", fill = NULL) +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out_dir, "mc_29_diagnostics.png"), gg, width = 8.5, height = 5, dpi = 120)
cat("\nGraphique -> mc_29_diagnostics.png\n")
