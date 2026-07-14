# Tests — Module 35 (optimiseurs modernes + regularisation NN). Ref : numDeriv.

test_that("Adam, RMSprop, momentum convergent vers l'optimum (quadratique mal conditionnee)", {
  b <- c(1, -2, 3); d <- c(1, 20, 0.5)
  g <- function(x) d * (x - b)
  expect_lt(max(abs(optim_adam(g, c(0, 0, 0), lr = 0.2, max_iter = 4000)$par - b)), 1e-3)
  expect_lt(max(abs(optim_rmsprop(g, c(0, 0, 0), lr = 0.1, max_iter = 6000)$par - b)), 1e-3)
  expect_lt(max(abs(optim_momentum(g, c(0, 0, 0), step = 0.02, max_iter = 6000)$par - b)), 1e-3)
})

test_that("Adam converge plus vite que la descente de gradient nue", {
  b <- c(1, -2, 3); d <- c(1, 20, 0.5); g <- function(x) d * (x - b)
  it_adam <- optim_adam(g, c(0, 0, 0), lr = 0.2, tol = 1e-6)$iters
  # GD au plus grand pas stable (2/L, L=20) : ~0.1
  gd <- function() { x <- c(0, 0, 0); for (t in 1:100000) { s <- 0.09 * g(x); x <- x - s
    if (max(abs(s)) < 1e-6) return(t) }; Inf }
  expect_lt(it_adam, gd())
})

test_that("dropout : identite en evaluation, esperance preservee en apprentissage", {
  set.seed(1); x <- rnorm(50000)
  expect_true(all(dropout(x, 0.5, training = FALSE)$out == x))
  expect_lt(abs(mean(dropout(x, 0.5, training = TRUE)$out) - mean(x)), 0.02)
})

test_that("batch_norm : sortie centree-reduite ; backward = numDeriv", {
  skip_if_not_installed("numDeriv")
  set.seed(2); X <- matrix(rnorm(30 * 4), 30, 4)
  bn <- batch_norm(X)
  expect_lt(max(abs(colMeans(bn$out))), 1e-10)             # moyennes ~ 0
  L <- function(v) sum(batch_norm(matrix(v, 30, 4))$out^2) / 2
  dX_home <- batch_norm_backward(bn$out, bn$cache)$dX      # dL/dout = out
  dX_num <- matrix(numDeriv::grad(L, as.numeric(X)), 30, 4)
  expect_lt(max(abs(dX_home - dX_num)), 1e-6)
})

test_that("batch_norm : gradients gamma/beta = numDeriv", {
  skip_if_not_installed("numDeriv")
  set.seed(3); X <- matrix(rnorm(20 * 3), 20, 3); g0 <- c(1.2, 0.8, 1.5); b0 <- c(0.1, -0.2, 0.3)
  bn <- batch_norm(X, g0, b0)
  bw <- batch_norm_backward(bn$out, bn$cache)
  Lg <- function(gv) sum(batch_norm(X, gv, b0)$out^2) / 2
  Lb <- function(bv) sum(batch_norm(X, g0, bv)$out^2) / 2
  expect_lt(max(abs(bw$dgamma - numDeriv::grad(Lg, g0))), 1e-6)
  expect_lt(max(abs(bw$dbeta - numDeriv::grad(Lb, b0))), 1e-6)
})
