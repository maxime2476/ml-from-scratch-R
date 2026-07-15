# Tests de conformité — Module 4 (régularisation). Références : forme fermée,
# MASS::lm.ridge, glmnet. Tolérances : 1e-8 (ridge exact), 1e-6 (lasso/glmnet).

skip_if_not_installed("glmnet")
skip_if_not_installed("MASS")

make_reg <- function(n = 150, p = 10, seed = 3, collinear = FALSE) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  if (collinear) X[, 2] <- X[, 1] + 0.05 * rnorm(n)
  beta0 <- c(2, 0, -1.5, 0, 0, 1, 0, 0, -0.7, 0)[seq_len(p)]
  y <- as.numeric(X %*% beta0 + rnorm(n))
  list(X = X, y = y, beta0 = beta0)
}

test_that("ridge_fit = forme fermée (X'X + lambda I)^-1 X'y (données brutes)", {
  d <- make_reg()
  for (lam in c(0.5, 5, 50)) {
    rf <- ridge_fit(d$X, d$y, lam, standardize = FALSE, intercept = FALSE)
    cf <- as.numeric(solve(crossprod(d$X) + lam * diag(ncol(d$X)), crossprod(d$X, d$y)))
    expect_equal(as.numeric(rf$beta), cf, tolerance = 1e-8)
  }
})

test_that("ridge_fit reproduit MASS::lm.ridge (intercept + pentes)", {
  d <- make_reg()
  for (lam in c(1, 10, 100)) {
    rf <- ridge_fit(d$X, d$y, lam, standardize = TRUE, intercept = TRUE)
    lr <- MASS::lm.ridge(d$y ~ d$X, lambda = lam)
    expect_equal(as.numeric(rf$coefficients), as.numeric(coef(lr)), tolerance = 1e-8)
  }
})

test_that("ridge rétrécit vers 0 quand lambda croît", {
  d <- make_reg()
  n0 <- sqrt(sum(ridge_fit(d$X, d$y, 1)$beta^2))
  n1 <- sqrt(sum(ridge_fit(d$X, d$y, 1e4)$beta^2))
  expect_lt(n1, n0)
})

test_that("soft_threshold : opérateur de seuillage doux (éq. 4.9)", {
  expect_equal(soft_threshold(c(-3, -0.5, 0, 0.5, 3), 1), c(-2, 0, 0, 0, 2))
  expect_equal(soft_threshold(5, 0), 5)                 # lambda=0 : identité
})

test_that("lasso_fit reproduit glmnet (alpha=1) sur design bien conditionné", {
  d <- make_reg()
  st <- mlfromscratch:::.standardize(d$X); Xs <- st$Xs; yc <- d$y - mean(d$y)
  n <- nrow(Xs)
  old <- glmnet::glmnet.control(); on.exit(do.call(glmnet::glmnet.control, old))
  glmnet::glmnet.control(thresh = 1e-14)
  for (lg in c(0.05, 0.2, 0.5)) {
    lf <- lasso_fit(Xs, yc, lambda = n * lg, standardize = FALSE, intercept = FALSE, tol = 1e-13)
    g1 <- glmnet::glmnet(Xs, yc, alpha = 1, lambda = lg, standardize = FALSE, intercept = FALSE)
    bg <- as.numeric(coef(g1))[-1]
    expect_equal(as.numeric(lf$beta), bg, tolerance = 1e-6)
    expect_identical(which(lf$beta != 0), which(bg != 0))   # même support
  }
})

test_that("lasso produit de la sparsité ; lambda=0 -> proche OLS", {
  d <- make_reg()
  # sparsité à lambda modéré
  lf <- lasso_fit(d$X, d$y, lambda = 30)
  expect_true(sum(lf$beta != 0) < ncol(d$X))
  # lambda ~ 0 : proche de l'OLS
  ols <- as.numeric(coef(lm(d$y ~ d$X)))
  lf0 <- lasso_fit(d$X, d$y, lambda = 1e-6, tol = 1e-12)
  expect_equal(as.numeric(lf0$coefficients), ols, tolerance = 1e-3)
})

test_that("ridge_bias_var : EQM analytique = EQM Monte Carlo", {
  set.seed(11)
  n <- 60; p <- 5
  X <- mlfromscratch:::.standardize(matrix(rnorm(n * p), n, p))$Xs
  beta_true <- c(2, -1, 0.5, 0, 1.5)
  sigma2 <- 1.5; lam <- 8
  bv <- ridge_bias_var(X, beta_true, sigma2, lam)
  # Monte Carlo : E[beta_hat] et Var, puis EQM = biais^2 + variance
  R <- 20000
  B <- matrix(0, R, p)
  for (r in seq_len(R)) {
    y <- as.numeric(X %*% beta_true + rnorm(n, sd = sqrt(sigma2)))
    B[r, ] <- ridge_fit(X, y, lam, standardize = FALSE, intercept = FALSE)$beta
  }
  bias_mc <- colMeans(B) - beta_true
  var_mc  <- apply(B, 2, var)
  mse_mc  <- sum(bias_mc^2) + sum(var_mc)
  expect_equal(bv$mse, mse_mc, tolerance = 0.05 * bv$mse)   # ~5 % (bruit MC)
})

test_that("ridge_bias_var : biais croît, variance décroît avec lambda", {
  set.seed(1)
  n <- 60; p <- 5
  X <- mlfromscratch:::.standardize(matrix(rnorm(n * p), n, p))$Xs
  bt <- c(2, -1, 0.5, 0, 1.5)
  bv_small <- ridge_bias_var(X, bt, 1, 1)
  bv_big   <- ridge_bias_var(X, bt, 1, 100)
  expect_gt(bv_big$bias2, bv_small$bias2)
  expect_lt(bv_big$variance, bv_small$variance)
})
