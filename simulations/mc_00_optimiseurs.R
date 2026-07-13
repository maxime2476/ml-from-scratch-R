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
trace_nesterov <- function(x0, step, K) {   # gradient accéléré + restart (éq. 0.22)
  x <- x0; y <- x0; lam <- 1; gap <- numeric(K)
  for (k in seq_len(K)) {
    g <- grad(y); x_new <- y - step * g
    if (sum(g * (x_new - x)) > 0) lam <- 1
    lam_new <- (1 + sqrt(1 + 4 * lam^2)) / 2
    y <- x_new + ((lam - 1) / lam_new) * (x_new - x)
    x <- x_new; lam <- lam_new; gap[k] <- fval(x) - fstar
  }
  data.frame(iter = seq_len(K), gap = pmax(gap, 1e-16), methode = "Nesterov (accéléré)")
}
trace_lbfgs <- function(x0, K) {            # L-BFGS instrumenté (mémoire m=10)
  x <- x0; g <- grad(x); S <- list(); Y <- list(); gap <- numeric(K)
  for (k in seq_len(K)) {
    q <- g; ns <- length(S); al <- numeric(ns); rho <- numeric(ns)
    for (i in rev(seq_len(ns))) { rho[i] <- 1/sum(Y[[i]]*S[[i]]); al[i] <- rho[i]*sum(S[[i]]*q); q <- q - al[i]*Y[[i]] }
    gam <- if (ns > 0) sum(S[[ns]]*Y[[ns]])/sum(Y[[ns]]*Y[[ns]]) else 1; r <- gam*q
    for (i in seq_len(ns)) { be <- rho[i]*sum(Y[[i]]*r); r <- r + S[[i]]*(al[i]-be) }
    d <- -r; slope <- sum(g*d); a <- 1; fx <- fval(x)
    while (fval(x + a*d) > fx + 1e-4*a*slope && a > 1e-12) a <- a/2
    xn <- x + a*d; gn <- grad(xn); s <- xn-x; yv <- gn-g
    if (sum(yv*s) > 1e-10) { S <- c(S,list(s)); Y <- c(Y,list(yv)); if(length(S)>10){S<-S[-1];Y<-Y[-1]} }
    x <- xn; g <- gn; gap[k] <- fval(x) - fstar
  }
  data.frame(iter = seq_len(K), gap = pmax(gap, 1e-16), methode = "L-BFGS (quasi-Newton)")
}

x0 <- rep(5, d)
K <- 60
trA <- rbind(trace_gd(x0, 1 / Lc, K), trace_nesterov(x0, 1 / Lc, K),
             trace_lbfgs(x0, K), trace_newton(x0, K))

cat(sprintf("=== (A) Après K=%d itérations : écart f_K - f* et it. pour < 1e-10 ===\n", K))
for (m in unique(trA$methode)) {
  sub <- trA[trA$methode == m, ]
  hit <- which(sub$gap < 1e-10)[1]
  cat(sprintf("  %-22s : gap %.2e  (< 1e-10 en %s it.)\n",
              m, tail(sub$gap, 1), if (is.na(hit)) paste0(">", K) else hit))
}
cat("Conditionnement kappa(A) =", round(Lc / min(eigen(A, only.values = TRUE)$values), 1), "\n")

pA <- ggplot(trA, aes(iter, gap, colour = methode)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  scale_y_log10() +
  labs(title = "Convergence : gradient vs Nesterov vs L-BFGS vs Newton",
       subtitle = "Gradient O(1/k) ; Nesterov accéléré ; L-BFGS quasi-Newton ; Newton quadratique",
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
