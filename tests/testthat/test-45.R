# Tests — Module 45 (GARCH). Reference : tseries::garch + recuperation des vrais parametres.

sim_garch <- function(n = 2000, w = 0.05, al = 0.1, be = 0.85, seed = 1) {
  set.seed(seed); x <- numeric(n); s2 <- numeric(n); s2[1] <- w / (1 - al - be)
  for (t in 2:n) { s2[t] <- w + al * x[t - 1]^2 + be * s2[t - 1]; x[t] <- rnorm(1) * sqrt(s2[t]) }
  x
}

test_that("GARCH(1,1) : parametres ~ tseries::garch", {
  skip_if_not_installed("tseries")
  x <- sim_garch()
  m <- garch_fit(x)
  gr <- tseries::garch(x, order = c(1, 1), trace = FALSE)
  expect_lt(abs(m$omega - coef(gr)[1]), 0.02)
  expect_lt(abs(m$alpha - coef(gr)[2]), 0.02)
  expect_lt(abs(m$beta - coef(gr)[3]), 0.02)
})

test_that("GARCH recupere les vrais parametres", {
  x <- sim_garch(n = 3000, w = 0.05, al = 0.1, be = 0.85)
  m <- garch_fit(x)
  expect_lt(abs(m$alpha - 0.1), 0.05)
  expect_lt(abs(m$beta - 0.85), 0.1)
  expect_lt(m$persistence, 1)                              # stationnarite (alpha+beta<1)
})

test_that("ARCH-LM detecte le regroupement de volatilite", {
  x_garch <- sim_garch()                                   # avec effet ARCH
  set.seed(2); x_iid <- rnorm(2000)                        # bruit blanc (pas d'ARCH)
  expect_lt(arch_lm_test(x_garch, 5)$p_value, 0.01)        # detecte
  expect_gt(arch_lm_test(x_iid, 5)$p_value, 0.05)          # ne detecte pas (correct)
})

test_that("la volatilite conditionnelle suit la volatilite realisee", {
  x <- sim_garch(n = 2000); m <- garch_fit(x)
  # correlation entre sigma_t^2 estime et |x_t| lisse (volatilite locale)
  local_vol <- stats::filter(x^2, rep(1 / 20, 20), sides = 1)
  ok <- !is.na(local_vol)
  expect_gt(cor(m$sigma[ok]^2, local_vol[ok]), 0.5)
})
