# =============================================================================
# Monte Carlo — Module 7 : KNN et fléau de la dimension
#  (A) Biais-variance du KNN selon k : variance ~ sigma²/k, biais croissant ;
#      l'EQM de test est en U (éq. 7.2).
#  (B) Fléau de la dimension : (i) concentration des distances (contraste et
#      cv(d²) -> 0 selon p, éq. 7.4-7.5) ; (ii) dégradation de la précision du
#      KNN quand on ajoute des dimensions de bruit.
# =============================================================================

source("R/07_knn.R")
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) Biais-variance du KNN (régression 1D)
# =============================================================================
f <- function(x) sin(3 * x)                     # vraie fonction
sigma <- 0.4; n <- 120
x_test <- seq(0.1, 0.9, length.out = 40)
f_test <- f(x_test)
ks <- c(1, 2, 3, 5, 8, 12, 20, 35, 60)
R <- 500

bvA <- data.frame(k = ks, biais2 = NA, variance = NA, eqm = NA)
for (ik in seq_along(ks)) {
  k <- ks[ik]
  preds <- matrix(NA_real_, R, length(x_test))
  for (r in seq_len(R)) {
    x <- runif(n); y <- f(x) + rnorm(n, sd = sigma)
    preds[r, ] <- knn_regression(matrix(x, n, 1), y, matrix(x_test, ncol = 1), k)
  }
  fbar <- colMeans(preds)
  bvA$biais2[ik]   <- mean((fbar - f_test)^2)
  bvA$variance[ik] <- mean(apply(preds, 2, var))
  bvA$eqm[ik]      <- bvA$biais2[ik] + bvA$variance[ik] + sigma^2
}
cat("=== (A) Biais-variance du KNN selon k ===\n")
print(round(bvA, 4), row.names = FALSE)
k_opt <- ks[which.min(bvA$eqm)]
cat(sprintf("k* (EQM min) = %d ; variance ~ sigma²/k (sigma²=%.2f)\n", k_opt, sigma^2))

longA <- reshape(bvA, varying = c("biais2", "variance", "eqm"), v.names = "valeur",
                 timevar = "composante", times = c("biais²", "variance", "EQM"),
                 direction = "long")
ggA <- ggplot(longA, aes(k, valeur, colour = composante)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  geom_vline(xintercept = k_opt, linetype = "dashed", colour = "grey40") +
  scale_x_log10() +
  labs(title = "KNN : décomposition biais-variance selon k",
       subtitle = "Petit k : variance élevée ; grand k : biais élevé ; EQM en U",
       x = "k (échelle log)", y = "valeur", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_07_biais_variance.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B-i) Concentration des distances selon p
# =============================================================================
ps <- c(2, 5, 10, 20, 50, 100, 200, 500)
concB <- data.frame(p = ps, contrast = NA, cv_d2 = NA)
for (i in seq_along(ps)) {
  dc <- distance_concentration(400, ps[i], seed = 1)
  concB$contrast[i] <- dc$contrast
  concB$cv_d2[i] <- dc$cv_d2
}
cat("\n=== (B-i) Concentration des distances selon p ===\n")
print(round(concB, 4), row.names = FALSE)
cat("Le contraste (Dmax-Dmin)/Dmin et cv(d²) tendent vers 0 (~1/sqrt(p)).\n")

ggBi <- ggplot(concB, aes(p)) +
  geom_line(aes(y = contrast, colour = "(Dmax−Dmin)/Dmin"), linewidth = 1) +
  geom_point(aes(y = contrast, colour = "(Dmax−Dmin)/Dmin")) +
  geom_line(aes(y = cv_d2, colour = "cv(d²) ~ 1/sqrt(p)"), linewidth = 1) +
  geom_point(aes(y = cv_d2, colour = "cv(d²) ~ 1/sqrt(p)")) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Fléau de la dimension : concentration des distances",
       subtitle = "Les voisins deviennent équidistants : le contraste s'effondre quand p croît",
       x = "dimension p (log)", y = "valeur (log)", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_07_concentration.png"), ggBi, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B-ii) Dégradation de la précision du KNN avec des dimensions de bruit
# =============================================================================
# Seules les 2 premières dimensions portent le signal ; on ajoute du bruit.
ps2 <- c(2, 5, 10, 25, 50, 100)
K <- 7; ntr <- 300; nte <- 300; Rb <- 40
accB <- data.frame(p = ps2, accuracy = NA)
for (i in seq_along(ps2)) {
  p <- ps2[i]; acc <- 0
  for (r in seq_len(Rb)) {
    Xtr <- matrix(rnorm(ntr * p), ntr, p); Xte <- matrix(rnorm(nte * p), nte, p)
    lab <- function(X) factor(ifelse(X[, 1] + X[, 2] > 0, "pos", "neg"))
    ytr <- lab(Xtr); yte <- lab(Xte)             # dépend SEULEMENT de x1, x2
    pred <- knn_classify(Xtr, ytr, Xte, K)
    acc <- acc + mean(pred == yte)
  }
  accB$accuracy[i] <- acc / Rb
}
cat("\n=== (B-ii) Précision KNN (k=7) avec dimensions de bruit ajoutées ===\n")
print(round(accB, 3), row.names = FALSE)
cat("Le signal est dans x1,x2 ; ajouter du bruit noie les distances et dégrade le KNN.\n")

ggBii <- ggplot(accB, aes(p, accuracy)) +
  geom_line(linewidth = 1, colour = "#2166ac") + geom_point(size = 2, colour = "#2166ac") +
  scale_x_log10() +
  labs(title = "Fléau de la dimension : le KNN se dégrade avec les variables de bruit",
       subtitle = "Signal dans x1, x2 seulement ; les dimensions inutiles diluent le voisinage",
       x = "dimension totale p (log)", y = "précision de test (k=7)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_07_knn_degradation.png"), ggBii, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir,
    "/mc_07_biais_variance.png, mc_07_concentration.png, mc_07_knn_degradation.png\n", sep = "")
