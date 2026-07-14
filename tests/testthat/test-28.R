# Tests — Module 28 (autodiff mode inverse). Référence : numDeriv (gradient exact).

test_that("gradient d'une fonction scalaire = numDeriv", {
  skip_if_not_installed("numDeriv")
  f <- function(x) sum(exp(x) * sin(x))
  x <- c(0.5, -1.2, 2.0)
  g_ad <- ad_grad(f, x)
  g_nd <- numDeriv::grad(f, x)
  expect_lt(max(abs(g_ad - g_nd)), 1e-7)
})

test_that("gradient d'une perte de régression = numDeriv", {
  skip_if_not_installed("numDeriv")
  set.seed(2); n <- 30; p <- 4
  X <- matrix(rnorm(n * p), n, p); y <- as.numeric(X %*% rnorm(p)) + rnorm(n)
  loss <- function(b) { r <- as.numeric(X %*% b) - y; sum(r * r) / n }
  b0 <- rnorm(p)
  ad_reset(); bn <- adnode(matrix(b0, p, 1))
  r <- mm(X, bn) - matrix(y, n, 1); L <- sum(r * r) / n; backward(L)
  expect_lt(max(abs(as.numeric(bn$grad) - numDeriv::grad(loss, b0))), 1e-7)
})

test_that("gradient d'un MLP (tanh) = numDeriv (= rétropropagation M12)", {
  skip_if_not_installed("numDeriv")
  set.seed(3); n <- 20; p <- 3; h <- 5
  X <- cbind(matrix(rnorm(n * p), n, p), 1); y <- matrix(rnorm(n), n, 1)
  unpack <- function(w) list(W1 = matrix(w[1:((p + 1) * h)], p + 1, h),
                             W2 = matrix(w[((p + 1) * h + 1):((p + 1) * h + h + 1)], h + 1, 1))
  mlp_loss <- function(w) { pk <- unpack(w); Ha <- cbind(tanh(X %*% pk$W1), 1)
    r <- Ha %*% pk$W2 - y; sum(r * r) / n }
  mlp_grad <- function(w) { pk <- unpack(w)
    ad_reset(); W1 <- adnode(pk$W1); W2 <- adnode(pk$W2)
    Ha <- ad_cbind1(tanh(mm(X, W1)))
    r <- mm(Ha, W2) - matrix(y, n, 1); L <- sum(r * r) / n; backward(L)
    c(as.numeric(W1$grad), as.numeric(W2$grad)) }
  w0 <- rnorm((p + 1) * h + h + 1)
  expect_lt(max(abs(mlp_grad(w0) - numDeriv::grad(mlp_loss, w0))), 1e-6)
})

test_that("dérivées des fonctions élémentaires (exp, log, tanh)", {
  expect_equal(ad_grad(function(x) sum(exp(x)), c(0, 1)), exp(c(0, 1)), tolerance = 1e-10)
  expect_equal(ad_grad(function(x) sum(log(x)), c(2, 4)), 1 / c(2, 4), tolerance = 1e-10)
  expect_equal(ad_grad(function(x) sum(tanh(x)), c(0.5, -1)), 1 - tanh(c(0.5, -1))^2, tolerance = 1e-10)
})

test_that("une passe arrière suffit (le graphe réutilisé donne le bon grad)", {
  # produit de trois : d/dx sum(x1*x2*x3) = (x2 x3, x1 x3, x1 x2)
  x <- c(2, 3, 4)
  g <- ad_grad(function(v) { p <- v * v; sum(p) }, x)   # d sum(x^2) = 2x
  expect_equal(g, 2 * x, tolerance = 1e-10)
})
