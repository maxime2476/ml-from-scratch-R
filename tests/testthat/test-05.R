# Tests de conformité — Module 5 (IV / 2SLS). Référence : AER::ivreg.
# Tolérance 1e-8.

skip_if_not_installed("AER")

make_iv <- function(n = 500, seed = 4) {
  set.seed(seed)
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n)
  x1 <- rnorm(n)
  xend <- 0.8 * z1 + 0.6 * z2 + w + rnorm(n)
  u <- 0.5 * w + rnorm(n)
  y <- 1 + 2 * x1 - 1.5 * xend + u
  data.frame(y, x1, xend, z1, z2)
}

test_that("2SLS sur-identifié reproduit AER::ivreg (coef, vcov, se, sigma2)", {
  d <- make_iv()
  X <- cbind(1, d$x1, d$xend); colnames(X) <- c("(Intercept)", "x1", "xend")
  Z <- cbind(1, d$x1, d$z1, d$z2); colnames(Z) <- c("(Intercept)", "x1", "z1", "z2")
  fit <- tsls_fit(d$y, X, Z)
  iv  <- AER::ivreg(y ~ x1 + xend | x1 + z1 + z2, data = d)
  expect_equal(as.numeric(fit$coefficients), as.numeric(coef(iv)), tolerance = 1e-8)
  expect_equal(unname(fit$vcov), unname(vcov(iv)), tolerance = 1e-8)
  expect_equal(unname(fit$se), unname(sqrt(diag(vcov(iv)))), tolerance = 1e-8)
  expect_equal(fit$sigma2, summary(iv)$sigma^2, tolerance = 1e-8)
})

test_that("2SLS juste-identifié = (Z'X)^-1 Z'y et = ivreg", {
  d <- make_iv()
  X <- cbind(1, d$x1, d$xend)
  Z <- cbind(1, d$x1, d$z1)
  fit <- tsls_fit(d$y, X, Z)
  beta_iv <- as.numeric(solve(crossprod(Z, X), crossprod(Z, d$y)))  # éq. 5.3
  expect_equal(as.numeric(fit$coefficients), beta_iv, tolerance = 1e-8)
  iv <- AER::ivreg(y ~ x1 + xend | x1 + z1, data = d)
  expect_equal(as.numeric(fit$coefficients), as.numeric(coef(iv)), tolerance = 1e-8)
})

test_that("first_stage_F reproduit le F de comparaison de modèles (lm/anova)", {
  d <- make_iv()
  Z <- cbind(1, d$x1, d$z1, d$z2)
  fs <- first_stage_F(d$xend, Z, excluded = c(3, 4))
  a <- anova(lm(xend ~ x1, d), lm(xend ~ x1 + z1 + z2, d))
  expect_equal(fs$F, a$F[2], tolerance = 1e-8)
  expect_equal(fs$df1, 2L)
})

test_that("2SLS corrige le biais d'endogénéité (grand échantillon)", {
  d <- make_iv(n = 5000)
  X <- cbind(1, d$x1, d$xend)
  Z <- cbind(1, d$x1, d$z1, d$z2)
  b_2sls <- tsls_fit(d$y, X, Z)$coefficients[3]
  b_ols  <- coef(lm(y ~ x1 + xend, d))[3]
  # 2SLS proche de -1.5 ; OLS nettement biaisé (attiré vers 0 par le confondeur)
  expect_lt(abs(b_2sls - (-1.5)), 0.1)
  expect_gt(abs(b_ols - (-1.5)), abs(b_2sls - (-1.5)))
})

test_that("2SLS = OLS quand Z = X (pas d'endogénéité instrumentée)", {
  d <- make_iv()
  X <- cbind(1, d$x1, d$z1)          # tout exogène ici
  fit <- tsls_fit(d$y, X, X)          # instruments = régresseurs
  ols <- solve_ls_qr(X, d$y)$coefficients
  expect_equal(as.numeric(fit$coefficients), as.numeric(ols), tolerance = 1e-8)
})
