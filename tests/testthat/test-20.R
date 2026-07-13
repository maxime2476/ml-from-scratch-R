# Tests — Module 20 (régression quantile). Référence : quantreg::rq.

skip_if_not_installed("quantreg")

make_q <- function(n = 400, seed = 3) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(n) * (1 + 0.5 * d$x1^2)
  d
}

test_that("pinball_loss : valeurs et forme", {
  expect_equal(pinball_loss(2, 0.9), 2 * 0.9)          # u>0 -> tau*u
  expect_equal(pinball_loss(-2, 0.9), -2 * (0.9 - 1))  # u<0 -> (tau-1)u = 0.2
  expect_equal(pinball_loss(4, 0.5), 2)                 # 0.5*|u|
})

test_that("Prop. 20.1 : la régression intercept-seul = quantile empirique (par la perte)", {
  set.seed(1); y <- rgamma(300, 2, 1); d <- data.frame(y = y)
  for (tau in c(0.25, 0.5, 0.75)) {
    q <- qreg_fit(y ~ 1, d, tau = tau)$coefficients[1]
    # la perte pinball en q est <= celle au quantile empirique (q est le minimiseur)
    expect_lte(sum(pinball_loss(y - q, tau)),
               sum(pinball_loss(y - quantile(y, tau), tau)) + 1e-6)
  }
})

test_that("qreg_fit atteint la MÊME perte optimale que quantreg::rq", {
  d <- make_q()
  for (tau in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    mine <- qreg_fit(y ~ x1 + x2, d, tau = tau, maxit = 500)
    ref  <- quantreg::rq(y ~ x1 + x2, tau = tau, data = d)
    loss_ref <- sum(pinball_loss(residuals(ref), tau))
    expect_lt(mine$loss, loss_ref + 1e-3)              # même minimum (LP)
  }
})

test_that("qreg_fit = rq pour les coefficients (quantiles centraux)", {
  d <- make_q()
  for (tau in c(0.5, 0.75)) {
    mine <- qreg_fit(y ~ x1 + x2, d, tau = tau, maxit = 500)
    ref  <- quantreg::rq(y ~ x1 + x2, tau = tau, data = d)
    expect_equal(as.numeric(mine$coefficients), as.numeric(coef(ref)), tolerance = 1e-4)
  }
})

test_that("régression médiane (LAD) robuste aux valeurs aberrantes vs OLS", {
  d <- make_q(); d$y[1:20] <- d$y[1:20] + 50           # 20 outliers
  b_ols <- coef(lm(y ~ x1 + x2, d))["x1"]
  b_lad <- qreg_fit(y ~ x1 + x2, d, tau = 0.5)$coefficients["x1"]
  # LAD reste proche de 2 ; OLS est tiré loin
  expect_lt(abs(b_lad - 2), abs(b_ols - 2))
})
