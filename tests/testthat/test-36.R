# Tests — Module 36 (CNN). Reference : definition + numDeriv.

test_that("conv2d = definition (correlation croisee)", {
  set.seed(1); X <- matrix(rnorm(36), 6, 6); K <- array(rnorm(9), c(3, 3, 1)); b <- 0.5
  cv <- conv2d(X, K, b)
  expect_equal(cv$out[2, 3, 1], b + sum(X[2:4, 3:5] * K[, , 1]), tolerance = 1e-12)
  expect_equal(dim(cv$out), c(4, 4, 1))
})

test_that("conv2d_backward = numDeriv (dX, dK, db)", {
  skip_if_not_installed("numDeriv")
  set.seed(2); H <- 7; W <- 7; kh <- 3; kw <- 3; Fn <- 2
  X <- matrix(rnorm(H * W), H, W); K <- array(rnorm(kh * kw * Fn), c(kh, kw, Fn)); b <- rnorm(Fn)
  cv <- conv2d(X, K, b); bw <- conv2d_backward(cv$out, cv$cache)   # dL/dout = out
  LX <- function(v) sum(conv2d(matrix(v, H, W), K, b)$out^2) / 2
  LK <- function(v) sum(conv2d(X, array(v, c(kh, kw, Fn)), b)$out^2) / 2
  Lb <- function(v) sum(conv2d(X, K, v)$out^2) / 2
  expect_lt(max(abs(bw$dX - matrix(numDeriv::grad(LX, as.numeric(X)), H, W))), 1e-5)
  expect_lt(max(abs(bw$dK - array(numDeriv::grad(LK, as.numeric(K)), c(kh, kw, Fn)))), 1e-5)
  expect_lt(max(abs(bw$db - numDeriv::grad(Lb, b))), 1e-6)
})

test_that("max_pool2d = maximum par bloc ; backward = numDeriv", {
  skip_if_not_installed("numDeriv")
  P <- matrix(1:16, 4, 4, byrow = TRUE)
  expect_equal(max_pool2d(P, 2)$out, matrix(c(6, 8, 14, 16), 2, 2, byrow = TRUE))
  set.seed(3); X <- matrix(rnorm(16), 4, 4); mp <- max_pool2d(X, 2)
  dX <- max_pool2d_backward(mp$out, mp$cache)
  L <- function(v) sum(max_pool2d(matrix(v, 4, 4), 2)$out^2) / 2
  expect_lt(max(abs(dX - matrix(numDeriv::grad(L, as.numeric(X)), 4, 4))), 1e-6)
})

test_that("convolution : equivariance a la translation", {
  set.seed(4); X <- matrix(0, 10, 10); X[3:5, 3:5] <- 2       # motif localise
  K <- array(rnorm(9), c(3, 3, 1))
  y1 <- conv2d(X, K)$out[, , 1]
  Xs <- matrix(0, 10, 10); Xs[3:5, 5:7] <- 2                  # meme motif decale de 2 en colonnes
  y2 <- conv2d(Xs, K)$out[, , 1]
  # la reponse est identique, decalee de 2 colonnes
  expect_lt(max(abs(y1[, 1:6] - y2[, 3:8])), 1e-10)
})
