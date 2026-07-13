# Tests — Module 18 (GMM). Références : tsls_fit (M5), lm, package gmm.

make_iv18 <- function(n = 500, seed = 4) {
  set.seed(seed)
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n); x1 <- rnorm(n)
  xend <- 0.8 * z1 + 0.6 * z2 + w + rnorm(n)
  u <- 2 * w + rnorm(n) * exp(0.3 * z1)            # hétéroscédastique
  y <- 1 + 2 * x1 - 1.5 * xend + u
  list(y = y, x1 = x1, xend = xend, z1 = z1, z2 = z2,
       X = cbind(`(Intercept)` = 1, x1 = x1, xend = xend),
       Z = cbind(`(Intercept)` = 1, x1 = x1, z1 = z1, z2 = z2))
}

test_that("GMM linéaire 1-étape (W=(Z'Z)^-1) = 2SLS (Prop. 18.2)", {
  d <- make_iv18()
  g <- gmm_linear(d$y, d$X, d$Z, twostep = FALSE)
  tv <- tsls_fit(d$y, d$X, d$Z)
  expect_equal(as.numeric(g$coefficients), as.numeric(tv$coefficients), tolerance = 1e-8)
})

test_that("GMM juste-identifié (Z = X) = OLS", {
  d <- make_iv18()
  Zx <- cbind(1, d$x1, d$xend)
  g <- gmm_linear(d$y, d$X, Zx, twostep = FALSE)
  expect_equal(as.numeric(g$coefficients), as.numeric(coef(lm(d$y ~ d$x1 + d$xend))), tolerance = 1e-8)
})

test_that("GMM efficace à deux étapes ~ package gmm (coef, se, J)", {
  skip_if_not_installed("gmm")
  d <- make_iv18()
  dat <- data.frame(y = d$y, x1 = d$x1, xend = d$xend, z1 = d$z1, z2 = d$z2)
  g2 <- gmm_linear(d$y, d$X, d$Z, twostep = TRUE)
  gp <- gmm::gmm(y ~ x1 + xend, ~ x1 + z1 + z2, data = dat, wmatrix = "optimal")
  expect_lt(max(abs(g2$coefficients - coef(gp))), 0.05)
  expect_lt(max(abs(g2$se - sqrt(diag(vcov(gp))))), 0.02)
  expect_equal(g2$J_df, 1L)
  expect_lt(abs(g2$J - gmm::specTest(gp)$test[1]), 0.5)
})

test_that("gmm_fit générique récupère moyenne et variance (moments)", {
  set.seed(1); v <- rnorm(500, 3, 2)
  gmom <- function(th, dat) cbind(dat - th[1], (dat - th[1])^2 - th[2])
  gf <- gmm_fit(gmom, c(0, 1), v, twostep = TRUE)
  expect_equal(gf$coefficients[1], mean(v), tolerance = 0.05)
  expect_equal(gf$coefficients[2], mean((v - mean(v))^2), tolerance = 0.1)
  expect_equal(gf$J_df, 0L)                          # juste-identifié : m=k=2
})

test_that("sous hétéroscédasticité, la GMM efficace réduit la variance vs 2SLS", {
  d <- make_iv18(n = 800)
  g_2sls <- gmm_linear(d$y, d$X, d$Z, twostep = FALSE)   # var robuste du 2SLS
  g_eff  <- gmm_linear(d$y, d$X, d$Z, twostep = TRUE)    # GMM efficace
  # la GMM efficace a une variance <= 2SLS sur le coef endogène (à l'échantillon près)
  expect_lte(g_eff$se["xend"], g_2sls$se["xend"] * 1.02)
})
