# Tests — Module 29 (diagnostics). Références : lmtest, AER, tseries.

make_diag <- function(n = 200, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(n) * (1 + abs(d$x1))   # hétéroscédastique
  d
}

test_that("Breusch-Pagan et White = lmtest::bptest", {
  skip_if_not_installed("lmtest")
  d <- make_diag(); fit <- lm(y ~ x1 + x2, d)
  expect_equal(bp_test(y ~ x1 + x2, d)$statistic,
               as.numeric(lmtest::bptest(fit)$statistic), tolerance = 1e-8)
  expect_equal(white_test(y ~ x1 + x2, d)$statistic,
               as.numeric(lmtest::bptest(fit, ~ x1 + x2 + I(x1^2) + I(x2^2) + I(x1 * x2),
                                         data = d)$statistic),
               tolerance = 1e-8)
})

test_that("Durbin-Watson et Breusch-Godfrey = lmtest", {
  skip_if_not_installed("lmtest")
  d <- make_diag(); fit <- lm(y ~ x1 + x2, d)
  expect_equal(dw_test(y ~ x1 + x2, d)$statistic,
               as.numeric(lmtest::dwtest(fit)$statistic), tolerance = 1e-8)
  for (o in 1:3)
    expect_equal(bg_test(y ~ x1 + x2, d, order = o)$statistic,
                 as.numeric(lmtest::bgtest(fit, order = o)$statistic), tolerance = 1e-8)
})

test_that("RESET = lmtest::resettest et Jarque-Bera = tseries", {
  skip_if_not_installed("lmtest")
  d <- make_diag(); fit <- lm(y ~ x1 + x2, d)
  expect_equal(reset_test(y ~ x1 + x2, d)$statistic,
               as.numeric(lmtest::resettest(fit, power = 2:3, type = "fitted")$statistic),
               tolerance = 1e-8)
  skip_if_not_installed("tseries")
  expect_equal(jarque_bera(resid(fit))$statistic,
               as.numeric(tseries::jarque.bera.test(resid(fit))$statistic), tolerance = 1e-8)
})

test_that("Durbin-Wu-Hausman et Sargan = AER::ivreg diagnostics", {
  skip_if_not_installed("AER")
  set.seed(7); n <- 300
  x1 <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n)
  xe <- 0.6 * z1 + 0.5 * z2 + 0.4 * x1 + rnorm(n)
  yy <- 1 + 2 * xe + x1 + rnorm(n)
  X <- cbind(1, xe, x1); Z <- cbind(1, z1, z2, x1)
  dg <- summary(AER::ivreg(yy ~ xe + x1 | z1 + z2 + x1), diagnostics = TRUE)$diagnostics
  expect_equal(dwh_test(yy, X, Z, endog = 2)$statistic,
               as.numeric(dg["Wu-Hausman", "statistic"]), tolerance = 1e-6)
  expect_equal(sargan_test(yy, X, Z)$statistic,
               as.numeric(dg["Sargan", "statistic"]), tolerance = 1e-6)
})

test_that("FGLS réduit l'hétéroscédasticité et reste sans biais", {
  d <- make_diag(n = 500)
  fg <- fgls(y ~ x1 + x2, d)
  expect_lt(max(abs(fg$coefficients - c(1, 2, -1))), 0.25)   # ~ vraies valeurs
  X <- model.matrix(lm(y ~ x1 + x2, d))
  rw <- as.numeric((d$y - X %*% fg$coefficients) * sqrt(fg$weights))
  bp_w <- 500 * summary(lm(I(rw^2) ~ d$x1 + d$x2))$r.squared
  expect_lt(bp_w, bp_test(y ~ x1 + x2, d)$statistic)         # BP pondéré < BP OLS
})
