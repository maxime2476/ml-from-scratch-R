# Tests — Module 14 (M-estimation, bayésien). Vérifications exactes reliant la
# théorie aux modules antérieurs.

skip_if_not_installed("sandwich")

test_that("sandwich M-estimation = HC0 de White (Th. 14.1, cas OLS)", {
  set.seed(3); n <- 200
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(n, sd = exp(0.4 * d$x1))
  fit <- ols_fit(y ~ x1 + x2, d)
  X <- fit$model_matrix
  V_m <- m_estimation_vcov(ols_psi(X, fit$residuals), crossprod(X) / n)
  expect_equal(unname(V_m), unname(vcov_hc(fit, "HC0")), tolerance = 1e-8)
})

test_that("postérieure du ridge : moyenne = estimateur ridge, cov = sigma2 (X'X+lI)^-1", {
  set.seed(5); n <- 150; p <- 4
  X <- matrix(rnorm(n * p), n, p); beta <- c(1, -2, 0.5, 0)
  y <- as.numeric(X %*% beta + rnorm(n))
  Xc <- scale(X, scale = FALSE); yc <- y - mean(y)
  lam <- 7; s2 <- 1.3
  post <- ridge_posterior(Xc, yc, lam, s2)
  beta_ridge <- as.numeric(solve(crossprod(Xc) + lam * diag(p), crossprod(Xc, yc)))
  expect_equal(post$mean, beta_ridge, tolerance = 1e-8)                 # Prop. 14.2
  expect_equal(post$cov, s2 * solve(crossprod(Xc) + lam * diag(p)), tolerance = 1e-8)
})

test_that("sous homoscédasticité, le sandwich se réduit à la variance classique", {
  set.seed(1); n <- 400
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + d$x1 - 0.5 * d$x2 + rnorm(n)          # erreurs homoscédastiques
  fit <- ols_fit(y ~ x1 + x2, d)
  X <- fit$model_matrix
  V_m <- m_estimation_vcov(ols_psi(X, fit$residuals), crossprod(X) / n)
  # proche de sigma^2 (X'X)^-1 (classique) : ratio des diagonales ~ 1
  ratio <- diag(V_m) / diag(fit$vcov)
  expect_true(all(abs(ratio - 1) < 0.15))
})
