# Tests de conformité — Module 0 (algèbre linéaire, optimiseurs).
# Références : qr()/qr.R, chol(), svd(), solve(), lm(). Tolérance 1e-8.

make_data <- function(n = 60, p = 5, seed = 123) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  beta <- rnorm(p)
  y <- as.numeric(X %*% beta + rnorm(n))
  list(X = X, y = y, beta = beta)
}

test_that("qr_householder reconstruit X et Q est orthogonale", {
  d <- make_data()
  dec <- qr_householder(d$X)
  expect_equal(dec$Q %*% dec$R, d$X, tolerance = 1e-8)
  expect_equal(crossprod(dec$Q), diag(nrow(d$X)), tolerance = 1e-8)
  # R triangulaire supérieure (sous-diagonale nulle)
  expect_true(all(abs(dec$R[lower.tri(dec$R)]) < 1e-8))
})

test_that("solve_ls_qr reproduit lm() : coefficients, résidus, RSS", {
  d <- make_data()
  fit <- lm(d$y ~ d$X - 1)
  res <- solve_ls_qr(d$X, d$y)
  expect_equal(as.numeric(res$coefficients), as.numeric(coef(fit)),
               tolerance = 1e-8)
  expect_equal(res$fitted, as.numeric(fitted(fit)), tolerance = 1e-8)
  expect_equal(res$rss, sum(residuals(fit)^2), tolerance = 1e-8)
})

test_that("chol_crout reproduit chol() (au transposé près)", {
  d <- make_data()
  A <- crossprod(d$X)
  L <- chol_crout(A)
  expect_equal(L %*% t(L), A, tolerance = 1e-8)         # A = L L^T
  expect_equal(L, t(chol(A)), tolerance = 1e-8)          # chol() renvoie R = L^T
  expect_true(all(abs(L[upper.tri(L)]) < 1e-12))         # L triangulaire inf.
})

test_that("solve_ls_chol reproduit lm()", {
  d <- make_data()
  fit <- lm(d$y ~ d$X - 1)
  res <- solve_ls_chol(d$X, d$y)
  expect_equal(as.numeric(res$coefficients), as.numeric(coef(fit)),
               tolerance = 1e-8)
})

test_that("substitutions triangulaires résolvent les systèmes", {
  set.seed(7)
  n <- 8
  L <- matrix(0, n, n); L[lower.tri(L, diag = TRUE)] <- rnorm(n * (n + 1) / 2)
  diag(L) <- abs(diag(L)) + 1
  b <- rnorm(n)
  expect_equal(forward_substitution(L, b), as.numeric(solve(L, b)),
               tolerance = 1e-8)
  U <- t(L)
  expect_equal(back_substitution(U, b), as.numeric(solve(U, b)),
               tolerance = 1e-8)
})

test_that("svd_tools : valeurs singulières, rang, conditionnement, pseudo-inverse", {
  d <- make_data()
  s <- svd_tools(d$X)
  ref <- svd(d$X)
  expect_equal(s$d, ref$d, tolerance = 1e-8)
  expect_equal(s$rank, ncol(d$X))
  expect_equal(s$kappa, ref$d[1] / ref$d[length(ref$d)], tolerance = 1e-8)
  # Axiomes de Moore-Penrose : X X^+ X = X et X^+ X X^+ = X^+
  Xp <- s$pinv
  expect_equal(d$X %*% Xp %*% d$X, d$X, tolerance = 1e-8)
  expect_equal(Xp %*% d$X %*% Xp, Xp, tolerance = 1e-8)
  # plein rang : X^+ = (X'X)^{-1} X'
  expect_equal(Xp, solve(crossprod(d$X), t(d$X)), tolerance = 1e-8)
})

test_that("solve_ls_svd = lm() en plein rang ; norme minimale en rang déficient", {
  d <- make_data()
  fit <- lm(d$y ~ d$X - 1)
  res <- solve_ls_svd(d$X, d$y)
  expect_equal(as.numeric(res$coefficients), as.numeric(coef(fit)),
               tolerance = 1e-8)
  expect_equal(res$rank, ncol(d$X))

  # colonne dupliquée : rang p, ajustement identique, solution de norme minimale
  Xd <- cbind(d$X, d$X[, 2])
  rd <- solve_ls_svd(Xd, d$y)
  expect_equal(rd$rank, ncol(d$X))
  expect_equal(rd$rss, sum(residuals(fit)^2), tolerance = 1e-8)
  # norme minimale : <= toute autre solution atteignant le même ajustement.
  # Une solution "naïve" met tout le poids sur la 1re copie, 0 sur la 2nde.
  beta_naif <- c(coef(fit), 0)
  expect_lte(sqrt(sum(rd$coefficients^2)), sqrt(sum(beta_naif^2)) + 1e-8)
})

test_that("optim_gd et optim_newton minimisent une quadratique SPD", {
  set.seed(11)
  A <- crossprod(matrix(rnorm(16), 4)) + diag(4)
  b <- rnorm(4)
  xstar <- as.numeric(solve(A, b))
  grad <- function(x) as.numeric(A %*% x - b)
  hess <- function(x) A
  L <- max(eigen(A, only.values = TRUE)$values)
  g <- optim_gd(grad, rep(0, 4), step = 1 / L, max_iter = 1e5, tol = 1e-12)
  nt <- optim_newton(grad, hess, rep(0, 4))
  expect_equal(g$par, xstar, tolerance = 1e-6)
  expect_equal(nt$par, xstar, tolerance = 1e-8)
  expect_lte(nt$iter, 5)   # convergence quadratique : très peu d'itérations
})

test_that("optim_cd résout une quadratique via argmin coordonnée fermé", {
  set.seed(13)
  A <- crossprod(matrix(rnorm(16), 4)) + diag(4)
  b <- rnorm(4)
  xstar <- as.numeric(solve(A, b))
  argmin_coord <- function(x, j) (b[j] - sum(A[j, -j] * x[-j])) / A[j, j]
  cd <- optim_cd(argmin_coord, rep(0, 4), tol = 1e-12)
  expect_equal(cd$par, xstar, tolerance = 1e-8)
})

test_that("optim_sgd approche la solution OLS", {
  d <- make_data(n = 200, p = 3)
  b_lm <- as.numeric(coef(lm(d$y ~ d$X - 1)))
  grad_i <- function(beta, idx) {
    Xi <- d$X[idx, , drop = FALSE]
    as.numeric(crossprod(Xi, as.numeric(Xi %*% beta) - d$y[idx])) / length(idx)
  }
  sg <- optim_sgd(grad_i, rep(0, 3), n = 200, batch = 10,
                  step_fun = function(t) 0.1 / (1 + 0.01 * t),
                  epochs = 400, seed = 1)
  expect_equal(sg$par, b_lm, tolerance = 1e-2)  # stochastique : tolérance lâche
})

test_that("optim_nesterov et optim_lbfgs minimisent une quadratique SPD", {
  set.seed(21)
  A <- crossprod(matrix(rnorm(64), 8)) + diag(8)
  b <- rnorm(8); xstar <- as.numeric(solve(A, b))
  f <- function(x) 0.5 * sum(x * (A %*% x)) - sum(b * x)
  grad <- function(x) as.numeric(A %*% x - b)
  L <- max(eigen(A, only.values = TRUE)$values)
  ne <- optim_nesterov(grad, rep(0, 8), step = 1 / L, max_iter = 1e5, tol = 1e-12)
  lb <- optim_lbfgs(grad, rep(0, 8), f = f, m = 10, max_iter = 200, tol = 1e-10)
  expect_equal(ne$par, xstar, tolerance = 1e-6)
  expect_equal(lb$par, xstar, tolerance = 1e-6)
})

test_that("Nesterov et L-BFGS convergent en MOINS d'itérations que le gradient", {
  set.seed(22); d <- 25
  M <- matrix(rnorm(d * d), d); A <- crossprod(M) + diag(d)   # mal conditionné
  b <- rnorm(d); f <- function(x) 0.5 * sum(x * (A %*% x)) - sum(b * x)
  grad <- function(x) as.numeric(A %*% x - b)
  L <- max(eigen(A, only.values = TRUE)$values); tol <- 1e-9
  gd <- optim_gd(grad, rep(0, d), step = 1 / L, max_iter = 1e5, tol = tol)
  ne <- optim_nesterov(grad, rep(0, d), step = 1 / L, max_iter = 1e5, tol = tol)
  lb <- optim_lbfgs(grad, rep(0, d), f = f, m = 10, max_iter = 1000, tol = 1e-8)
  expect_lt(ne$iter, gd$iter)                        # accélération de Nesterov
  expect_lt(lb$iter, gd$iter)                        # quasi-Newton plus rapide
})
