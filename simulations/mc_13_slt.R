# =============================================================================
# Monte Carlo / illustrations numériques — Module 13 : théorie de l'apprentissage
#  (A) Hoeffding empirique : fréquence des déviations vs borne (13.2).
#  (B) Complexité de Rademacher : classe linéaire (vs borne 13.5) et souches
#      d'arbre ; écart train/test observé vs borne de généralisation (13.4).
#  (C) Pulvérisation : 3 points de R^2 pulvérisés par les hyperplans, 4 non.
# =============================================================================

for (f in c("00_linalg", "01_ols", "03_glm_irls", "13_slt"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages({library(ggplot2); library(quadprog)})

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) Hoeffding empirique
# =============================================================================
n <- 120; p <- 0.5; R <- 50000
eps_grid <- seq(0.05, 0.30, by = 0.025)
Rhat <- colMeans(matrix(rbinom(n * R, 1, p), n, R))
tabA <- data.frame(eps = eps_grid,
                   freq = sapply(eps_grid, function(e) mean(abs(Rhat - p) >= e)),
                   borne = hoeffding_bound(n, eps_grid))
cat("=== (A) Hoeffding : P(|Rhat - R| >= eps) vs borne (n =", n, ") ===\n")
print(round(tabA, 4), row.names = FALSE)
cat("La fréquence empirique reste SOUS la borne pour tout eps (borne valide, et lâche).\n")

ggA <- ggplot(tabA, aes(eps)) +
  geom_line(aes(y = borne, colour = "borne 2e^{-2n eps^2}"), linewidth = 1) +
  geom_point(aes(y = freq, colour = "fréquence empirique"), size = 2) +
  geom_line(aes(y = freq, colour = "fréquence empirique")) +
  scale_y_log10() +
  labs(title = "Inégalité de Hoeffding : empirique vs borne",
       subtitle = paste0("n = ", n, ", perte de Bernoulli ; la borne majore la fréquence des déviations"),
       x = expression(epsilon), y = "probabilité (log)", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_13_hoeffding.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Complexité de Rademacher : linéaire vs souches ; écart de généralisation
# =============================================================================
# Rademacher empirique de la classe des souches (stumps) axis-aligned, sorties +-1.
stump_rademacher <- function(X, n_draws = 1500L, seed = NULL) {
  X <- as.matrix(X); n <- nrow(X); d <- ncol(X)
  if (!is.null(seed)) set.seed(seed)
  ords <- lapply(seq_len(d), function(j) order(X[, j]))
  vals <- numeric(n_draws)
  for (b in seq_len(n_draws)) {
    sigma <- sample(c(-1, 1), n, replace = TRUE)
    best <- 0
    for (j in seq_len(d)) {
      s <- sigma[ords[[j]]]
      pref <- c(0, cumsum(s))           # prefix sums
      tot <- pref[n + 1]
      # sup_k | tot - 2*pref_k | (seuil entre positions k et k+1, deux polarités)
      best <- max(best, max(abs(tot - 2 * pref)))
    }
    vals[b] <- best
  }
  mean(vals) / n
}

ns <- c(50L, 100L, 200L, 400L)
d <- 5
tabB <- data.frame(n = ns, rad_lineaire = NA, borne_lineaire = NA, rad_souches = NA)
for (i in seq_along(ns)) {
  X <- matrix(rnorm(ns[i] * d), ns[i], d)
  tabB$rad_lineaire[i]   <- empirical_rademacher_linear(X, B = 1, n_draws = 2000, seed = 1)
  tabB$borne_lineaire[i] <- rademacher_linear_bound(X, B = 1)
  tabB$rad_souches[i]    <- stump_rademacher(X, n_draws = 1500, seed = 1)
}
cat("\n=== (B1) Complexité de Rademacher empirique selon n (d =", d, ") ===\n")
print(round(tabB, 4), row.names = FALSE)
cat("Linéaire : sous la borne B*rho/sqrt(n) ; souches : complexité comparable, ~1/sqrt(n).\n")

# Écart de généralisation d'un classifieur linéaire (ERM logistique) vs borne
gap_experiment <- function(n, d = 5, reps = 300, delta = 0.05) {
  gaps <- numeric(reps); rad <- numeric(reps)
  for (r in seq_len(reps)) {
    beta <- c(0.8, rep(c(1, -1), length.out = d))
    Xtr <- matrix(rnorm(n * d), n, d)
    ytr <- rbinom(n, 1, plogis(cbind(1, Xtr) %*% beta))
    Xte <- matrix(rnorm(2000 * d), 2000, d)
    yte <- rbinom(2000, 1, plogis(cbind(1, Xte) %*% beta))
    fit <- tryCatch(glm_irls(y ~ ., data.frame(y = ytr, Xtr), "binomial"),
                    warning = function(w) NULL, error = function(e) NULL)
    if (is.null(fit)) next
    pr_tr <- as.numeric(cbind(1, Xtr) %*% fit$coefficients) > 0
    pr_te <- as.numeric(cbind(1, Xte) %*% fit$coefficients) > 0
    err_tr <- mean(pr_tr != ytr); err_te <- mean(pr_te != yte)
    gaps[r] <- err_te - err_tr
    rad[r]  <- empirical_rademacher_linear(Xtr, B = 1, n_draws = 400)
  }
  bound <- 2 * mean(rad) + sqrt(log(1 / delta) / (2 * n))
  c(gap_moyen = mean(gaps), gap_max = quantile(gaps, 0.95), borne = bound)
}
cat("\n=== (B2) Écart train/test d'un classifieur linéaire vs borne (13.4) ===\n")
gapTab <- as.data.frame(t(sapply(c(50L, 100L, 200L), function(nn) gap_experiment(nn))))
gapTab$n <- c(50, 100, 200)
print(round(gapTab[, c("n", "gap_moyen", "gap_max.95%", "borne")], 3), row.names = FALSE)
cat("L'écart observé reste SOUS la borne — qui tient mais est très lâche (typique).\n")

longB <- reshape(tabB, varying = c("rad_lineaire", "borne_lineaire", "rad_souches"),
                 v.names = "valeur", timevar = "quantite",
                 times = c("Rademacher linéaire", "borne linéaire", "Rademacher souches"),
                 direction = "long")
ggB <- ggplot(longB, aes(n, valeur, colour = quantite)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  labs(title = "Complexité de Rademacher : linéaire vs souches d'arbre",
       subtitle = "Décroissance en ~1/sqrt(n) ; la borne B*rho/sqrt(n) majore le cas linéaire",
       x = "n", y = "complexité de Rademacher empirique", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_13_rademacher.png"), ggB, width = 8, height = 5, dpi = 120)

# =============================================================================
# (C) Pulvérisation : 3 points oui, 4 non
# =============================================================================
P3 <- matrix(c(0, 0, 1, 0, 0.5, 1), 3, 2, byrow = TRUE)
P4 <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), 4, 2, byrow = TRUE)
cat("\n=== (C) Pulvérisation par les demi-plans de R^2 ===\n")
cat("3 points (position générale) pulvérisés :", shatters_hyperplane(P3),
    "-> VC >= 3\n")
cat("4 points (carré)             pulvérisés :", shatters_hyperplane(P4),
    "-> VC < 4, donc VC = 3\n")

# Illustration : l'étiquetage XOR des 4 points n'est pas linéairement séparable
dfC <- data.frame(x = P4[, 1], y = P4[, 2], label = factor(c(1, -1, -1, 1)))
ggC <- ggplot(dfC, aes(x, y, colour = label, shape = label)) +
  geom_point(size = 6) +
  geom_abline(slope = -1, intercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_abline(slope = -1, intercept = 1.5, linetype = "dashed", colour = "grey50") +
  labs(title = "VC-dim des demi-plans = 3 : l'étiquetage XOR de 4 points est irréalisable",
       subtitle = "Aucune droite ne sépare la diagonale (+) de l'anti-diagonale (−)",
       x = expression(x[1]), y = expression(x[2])) +
  coord_equal() + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_13_shattering.png"), ggC, width = 6, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir,
    "/mc_13_hoeffding.png, mc_13_rademacher.png, mc_13_shattering.png\n", sep = "")
