# =============================================================================
# Monte Carlo — Module 26 : la courbe de double descente
# On moyenne la courbe de risque (train/test) sur plusieurs jeux, sans puis avec
# ridge. Sans régularisation : pic à D=n puis SECONDE descente. Avec ridge : pic
# supprimé (courbe monotone). Illustre que la double descente est un artefact de
# l'interpolation NON régularisée.
# =============================================================================

for (f in c("00_linalg", "04_regularisation", "26_double_descent"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n <- 80; p <- 5; R <- 30
f0 <- function(X) sin(2 * X[, 1]) + 0.5 * X[, 2]^2 - X[, 3]
Ds <- c(5, 15, 30, 45, 60, 70, 76, 80, 84, 90, 100, 130, 180, 260, 400, 600)

acc <- function(lambda) {
  M <- matrix(0, length(Ds), 2)
  for (r in seq_len(R)) {
    Xtr <- matrix(rnorm(n * p), n, p); ytr <- f0(Xtr) + 0.3 * rnorm(n)
    Xte <- matrix(rnorm(1000 * p), 1000, p); yte <- f0(Xte) + 0.3 * rnorm(1000)
    cur <- double_descent_curve(Xtr, ytr, Xte, yte, Ds, gamma = 0.5, seed = 100 + r, lambda = lambda)
    M <- M + as.matrix(cur[, c("train_mse", "test_mse")])
  }
  M / R
}
m0 <- acc(0); m1 <- acc(2)

cat("=== Double descente (n =", n, ", moyenne sur", R, "jeux) ===\n\n")
cat(sprintf("%6s %12s %12s\n", "D", "test (min-norm)", "test (ridge)"))
for (i in seq_along(Ds))
  cat(sprintf("%6d %12.3f %12.3f%s\n", Ds[i], m0[i, 2], m1[i, 2],
              if (Ds[i] == n) "   <- seuil D=n" else ""))
cat(sprintf("\n=> min-norm : pic au seuil (test=%.1f) puis descente jusqu'à %.3f (< régime\n",
            m0[which(Ds == n), 2], min(m0[Ds > n, 2])))
cat(sprintf("   sous-paramétré, min=%.3f). Ridge : pas de pic (max=%.3f, monotone-ish).\n",
            min(m0[Ds < n, 2]), max(m1[, 2])))

df <- rbind(
  data.frame(D = Ds, mse = m0[, 2], regime = "test — min-norm (lambda=0)"),
  data.frame(D = Ds, mse = m1[, 2], regime = "test — ridge (lambda=2)"),
  data.frame(D = Ds, mse = m0[, 1], regime = "apprentissage — min-norm"))
gg <- ggplot(df, aes(D, mse, colour = regime)) +
  geom_vline(xintercept = n, linetype = "dashed", colour = "grey60") +
  annotate("text", x = n, y = max(m0[, 2]), label = "D = n", hjust = -0.1, size = 3.5) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  scale_y_log10() +
  labs(title = "Double descente : le risque de test rechute après le seuil d'interpolation",
       subtitle = "min-norm : pic à D=n puis seconde descente ; le ridge supprime le pic",
       x = "nombre de caractéristiques D (complexité)", y = "EQM (échelle log)", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_26_double_descent.png"), gg, width = 9, height = 5.5, dpi = 120)
cat("\nGraphique -> mc_26_double_descent.png\n")
