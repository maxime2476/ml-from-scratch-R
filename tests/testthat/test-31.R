# Tests — Module 31 (series temporelles). References : stats, tseries.

make_ts <- function(n = 300, seed = 1) { set.seed(seed)
  as.numeric(arima.sim(list(ar = c(0.6, -0.2)), n)) }

test_that("ACF et PACF = stats::acf / stats::pacf", {
  x <- make_ts()
  expect_equal(acf_ts(x, 8), as.numeric(acf(x, 8, plot = FALSE)$acf), tolerance = 1e-8)
  expect_equal(pacf_ts(x, 8), as.numeric(pacf(x, 8, plot = FALSE)$acf), tolerance = 1e-8)
})

test_that("AR Yule-Walker = stats::ar", {
  x <- make_ts()
  yh <- ar_yw(x, 2); yr <- ar(x, aic = FALSE, order.max = 2, method = "yule-walker")
  expect_equal(yh$ar, as.numeric(yr$ar), tolerance = 1e-8)
})

test_that("Ljung-Box = stats::Box.test", {
  x <- make_ts()
  expect_equal(ljung_box(x, 10)$statistic,
               as.numeric(Box.test(x, 10, type = "Ljung-Box")$statistic), tolerance = 1e-8)
})

test_that("ADF = tseries::adf.test (meme statistique)", {
  skip_if_not_installed("tseries")
  x <- make_ts()
  expect_equal(adf_test(x)$statistic,
               as.numeric(tseries::adf.test(x)$statistic), tolerance = 1e-6)
})

test_that("ADF distingue stationnaire (rejet) de racine unitaire (non-rejet)", {
  set.seed(2)
  stat <- as.numeric(arima.sim(list(ar = 0.3), 300))     # stationnaire
  rw <- cumsum(rnorm(300))                                # marche aleatoire (racine unitaire)
  expect_lt(adf_test(stat)$statistic, -3)                # tres negatif -> rejet
  expect_gt(adf_test(rw)$statistic, -3)                  # proche de 0 -> non rejet
})

test_that("ARMA-CSS recupere les vrais parametres et ~ arima(CSS)", {
  set.seed(5)
  x <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.4), 1000))
  fit <- arma_css(x, 1, 1)
  expect_lt(abs(fit$ar - 0.5), 0.1)                       # recupere ar
  expect_lt(abs(fit$ma - 0.4), 0.1)                       # recupere ma
  ar_ref <- arima(x, order = c(1, 0, 1), method = "CSS")
  expect_lt(abs(fit$ar - coef(ar_ref)[1]), 0.03)          # ~ arima
})
