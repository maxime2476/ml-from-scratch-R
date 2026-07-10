# =============================================================================
# Monte Carlo — Module 0 : comportement des optimiseurs
# Deux questions :
#  (A) Vitesse de convergence : gradient (O(1/k), linéaire si fortement convexe)
#      vs Newton (quadratique) sur une quadratique fortement convexe.
#  (B) Compromis pas/variance du SGD : pas constant (plancher de variance) vs
#      pas décroissant (Robbins-Monro, éq. 0.21 -> convergence).
# =============================================================================

source("R/00_linalg.R")   # exécuter depuis la racine du projet
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) Gradient vs Newton — quadratique f(x) = 1/2 x' A x - b' x
# =============================================================================
d <- 20
M <- matrix(rnorm(d * d), d)
A <- crossprod(M) + diag(d)                 # SPD
b <- rnorm(d)
xstar <- as.numeric(solve(A, b))
fstar <- 0.5 * as.numeric(t(xstar) %*% A %*% xstar) - sum(b * xstar)
fval  <- function(x) 0.5 * as.numeric(t(x) %*% A %*% x) - sum(b * x)
grad  <- function(x) as.numeric(A %*% x - b)
Lc <- max(eigen(A, only.values = TRUE)$values)   # constante de lissité

# Traces d'itérés (on ré-instrumente à la main pour capturer f_k - f*).
trace_gd <- function(x0, step, K) {
  x <- x0; gap <- numeric(K)
  for (k in seq_len(K)) { x <- x - step * grad(x); gap[k] <- fval(x) - fstar }
  data.frame(iter = seq_len(K), gap = pmax(gap, 1e-16), methode = "Gradient (t=1/L)")
}
trace_newton <- function(x0, K) {
  x <- x0; gap <- numeric(K)
  for (k in seq_len(K)) { x <- x - solve(A, grad(x)); gap[k] <- fval(x) - fstar }
  data.frame(iter = seq_len(K), gap = pmax(gap, 1e-16), methode = "Newton")
}

x0 <- rep(5, d)
K <- 60
trA <- rbind(trace_gd(x0, 1 / Lc, K), trace_newton(x0, K))

cat("=== (A) Itérations pour atteindre f_k - f* < 1e-10 ===\n")
for (m in unique(trA$methode)) {
  sub <- trA[trA$methode == m, ]
  hit <- which(sub$gap < 1e-10)[1]
  cat(sprintf("  %-18s : %s itérations\n", m, if (is.na(hit)) ">K" else hit))
}
cat("Conditionnement kappa(A) =", round(Lc / min(eigen(A, only.values = TRUE)$values), 1), "\n")

pA <- ggplot(trA, aes(iter, gap, colour = methode)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  scale_y_log10() +
  labs(title = "Convergence : gradient vs Newton (quadratique fortement convexe)",
       subtitle = "Gradient linéaire (droite en échelle log) ; Newton quadratique (chute brutale)",
       x = "itération k", y = expression(f(x[k]) - f^"*"), colour = "Méthode") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_00_optim_convergence.png"), pA,
       width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) SGD : pas constant vs pas décroissant
# Perte empirique de régression : f(beta) = 1/(2n) sum (x_i' beta - y_i)^2.
# =============================================================================
n <- 400; p <- 5
Xd <- cbind(1, matrix(rnorm(n * (p - 1)), n))
beta0 <- c(1, -2, 0.5, 0, 3)
yv <- as.numeric(Xd %*% beta0 + rnorm(n))
beta_ols <- as.numeric(solve(crossprod(Xd), crossprod(Xd, yv)))

grad_i <- function(beta, idx) {
  Xi <- Xd[idx, , drop = FALSE]
  as.numeric(crossprod(Xi, as.numeric(Xi %*% beta) - yv[idx])) / length(idx)
}

# SGD instrumenté : distance à l'optimum OLS à chaque update.
run_sgd_trace <- function(step_fun, epochs, batch = 1, seed = 1) {
  set.seed(seed)
  beta <- rep(0, p); t <- 0L; dist <- numeric(0)
  for (e in seq_len(epochs)) {
    idx <- sample.int(n)
    for (start in seq(1L, n, by = batch)) {
      t <- t + 1L
      bi <- idx[start:min(start + batch - 1L, n)]
      beta <- beta - step_fun(t) * grad_i(beta, bi)
      dist <- c(dist, sqrt(sum((beta - beta_ols)^2)))
    }
  }
  dist
}

epochs <- 60; batch <- 1
runs <- list(
  "Constant t=0.05"       = function(t) 0.05,
  "Constant t=0.005"      = function(t) 0.005,
  "Décroissant t0/(1+ct)" = function(t) 0.1 / (1 + 0.02 * t)
)
traceB <- do.call(rbind, lapply(names(runs), function(nm) {
  dvec <- run_sgd_trace(runs[[nm]], epochs, batch, seed = 7)
  data.frame(update = seq_along(dvec), dist = dvec, schema = nm)
}))

cat("\n=== (B) Distance finale à l'optimum OLS (moyenne 200 derniers updates) ===\n")
for (nm in names(runs)) {
  sub <- traceB[traceB$schema == nm, ]
  cat(sprintf("  %-24s : %.4e\n", nm, mean(tail(sub$dist, 200))))
}
cat("Lecture : pas constant -> plancher de variance ~ t ; pas décroissant -> tend vers 0.\n")

pB <- ggplot(traceB, aes(update, dist, colour = schema)) +
  geom_line(alpha = 0.8) +
  scale_y_log10() +
  labs(title = "SGD : compromis pas / variance",
       subtitle = "Pas constant = descente rapide puis plancher de variance ; pas décroissant (Robbins-Monro) = convergence",
       x = "mises à jour", y = expression("||"*beta[t] - hat(beta)[OLS]*"||"),
       colour = "Schéma de pas") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_00_optim_sgd.png"), pB,
       width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_00_optim_convergence.png, mc_00_optim_sgd.png\n", sep = "")
