# =============================================================================
# Monte Carlo / illustrations — Module 15 : interprétabilité
#  (A) Piège du PDP sous corrélation (deux variables corrélées à 0.9) : le PDP
#      évalue le modèle en des points hors du support conjoint (extrapolation).
#  (B) ICE vs PDP sous interaction : le PDP (plat) masque l'hétérogénéité que
#      l'ICE révèle.
#  (C) Convergence de l'estimateur de Shapley par permutations vers l'exact.
# =============================================================================

for (f in c("00_linalg", "08_cart", "09_bagging_rf", "15_interpretabilite"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages({library(ggplot2); library(MASS)})

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) Piège du PDP sous corrélation
# =============================================================================
gen_corr <- function(n, rho) {
  X <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, rho, rho, 1), 2))
  data.frame(x1 = X[, 1], x2 = X[, 2])
}
# Fraction "hors support" : distance de Mahalanobis des points synthétiques du
# PDP (grille_x1 x x2 observés) sous la loi jointe des données.
oos_fraction <- function(dat, grid) {
  S <- cov(dat); m <- colMeans(dat)
  synth <- expand.grid(x1 = grid, x2 = dat$x2)
  d2 <- mahalanobis(synth, m, S)                  # mahalanobis attend la covariance
  mean(d2 > qchisq(0.95, df = 2))                 # au-delà de l'ellipse à 95 %
}

n <- 400
cat("=== (A) Piège du PDP : fraction de points synthétiques hors support ===\n")
for (rho in c(0, 0.5, 0.9)) {
  d <- gen_corr(n, rho)
  grid <- seq(min(d$x1), max(d$x1), length.out = 25)
  cat(sprintf("rho=%.1f : %.1f %% des points (grille x1) x (x2 obs) hors de l'ellipse 95%%\n",
              rho, 100 * oos_fraction(d, grid)))
}
cat("Plus la corrélation est forte, plus le PDP évalue le modèle hors des données.\n")

# Visualisation : données réelles (rho=0.9) + points synthétiques du PDP
d9 <- gen_corr(n, 0.9)
grid <- seq(min(d9$x1), max(d9$x1), length.out = 12)
synth <- expand.grid(x1 = grid, x2 = d9$x2)
ggA <- ggplot() +
  geom_point(data = synth, aes(x1, x2), colour = "#d73027", alpha = 0.06, size = 0.8) +
  geom_point(data = d9, aes(x1, x2), colour = "#1b7837", alpha = 0.5, size = 1) +
  labs(title = "Piège du PDP : les points synthétiques sortent du support (rho=0.9)",
       subtitle = "Vert = données réelles (nuage diagonal) ; rouge = grille du PDP (remplit le carré)",
       x = expression(x[1]), y = expression(x[2])) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_15_pdp_piege.png"), ggA, width = 7, height = 5, dpi = 120)

# =============================================================================
# (B) ICE vs PDP sous interaction
# =============================================================================
set.seed(3)
Xb <- data.frame(x1 = runif(200, -2, 2), x2 = sample(c(-1, 1), 200, replace = TRUE))
# interaction : l'effet de x1 change de SIGNE selon x2
pf_int <- function(D) D$x1 * D$x2
ic <- ice(pf_int, Xb, "x1", grid_size = 30)
pente_pdp <- coef(lm(ic$pdp ~ ic$grid))[2]
cat("\n=== (B) ICE vs PDP sous interaction (effet de x1 dépend de x2) ===\n")
cat(sprintf("pente du PDP (moyenne) = %.3f (proche de 0 : effet apparent nul)\n", pente_pdp))
cat("Les courbes ICE ont des pentes +1 et -1 selon x2 : le PDP les moyenne à ~0.\n")

df_ice <- data.frame(
  x1 = rep(ic$grid, each = nrow(Xb)),
  val = as.numeric(t(ic$ice)),
  id = rep(seq_len(nrow(Xb)), times = length(ic$grid)),
  x2 = rep(Xb$x2, times = length(ic$grid)))
ggB <- ggplot() +
  geom_line(data = df_ice, aes(x1, val, group = id, colour = factor(x2)), alpha = 0.25) +
  geom_line(data = data.frame(x1 = ic$grid, pdp = ic$pdp), aes(x1, pdp),
            colour = "black", linewidth = 1.3) +
  labs(title = "ICE vs PDP sous interaction",
       subtitle = "Courbes ICE de pentes opposées (couleur = x2) ; PDP (noir) ~ plat les masque",
       x = expression(x[1]), y = "prédiction", colour = "x2") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_15_ice_pdp.png"), ggB, width = 8, height = 5, dpi = 120)

# =============================================================================
# (C) Convergence de Shapley par permutations
# =============================================================================
set.seed(5)
Xc <- data.frame(x1 = rnorm(200), x2 = rnorm(200), x3 = rnorm(200),
                 x4 = rnorm(200), x5 = rnorm(200))
beta <- c(2, -1.5, 1, -0.5, 0.3)
pf_lin <- function(D) as.numeric(as.matrix(D[, 1:5]) %*% beta)
phi_exact <- shapley_exact(pf_lin, Xc[1, ], Xc)
ns <- c(50, 100, 200, 500, 1000, 2000, 5000)
errC <- data.frame(n_samples = ns, err = sapply(ns, function(ns_)
  max(abs(shapley_permutation(pf_lin, Xc[1, ], Xc, n_samples = ns_, seed = 1) - phi_exact))))
cat("\n=== (C) Convergence de Shapley (permutations) vers l'exact ===\n")
print(round(errC, 4), row.names = FALSE)
cat("L'erreur décroît en ~1/sqrt(n_samples) (estimateur Monte Carlo).\n")

ggC <- ggplot(errC, aes(n_samples, err)) +
  geom_line(linewidth = 1, colour = "#2166ac") + geom_point(size = 2, colour = "#2166ac") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Shapley : convergence de l'estimateur par permutations",
       subtitle = "Erreur max vs valeur exacte ~ 1/sqrt(nombre de tirages)",
       x = "nombre de tirages (log)", y = "erreur max (log)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_15_shapley_conv.png"), ggC, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir,
    "/mc_15_pdp_piege.png, mc_15_ice_pdp.png, mc_15_shapley_conv.png\n", sep = "")
