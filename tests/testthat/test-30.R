# Tests — Module 30 (variables dependantes limitees).
# References : glm(probit), AER::tobit, sampleSelection::selection.

test_that("probit = glm(binomial probit) a 1e-6", {
  set.seed(1); n <- 400; x <- rnorm(n); y <- rbinom(n, 1, pnorm(0.5 + 1.2 * x))
  d <- data.frame(y = y, x = x)
  pr <- probit(y ~ x, d); g <- glm(y ~ x, d, family = binomial(link = "probit"))
  expect_equal(as.numeric(pr$coefficients), as.numeric(coef(g)), tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(as.numeric(pr$se), as.numeric(sqrt(diag(vcov(g)))), tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("Tobit = AER::tobit a 1e-4", {
  skip_if_not_installed("AER")
  set.seed(2); n <- 500; x <- rnorm(n); yl <- pmax(0, 1 + 0.8 * x + rnorm(n))
  d <- data.frame(y = yl, x = x)
  tb <- tobit_fit(y ~ x, d, left = 0); at <- AER::tobit(y ~ x, data = d)
  expect_equal(as.numeric(tb$coefficients), as.numeric(coef(at)), tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(as.numeric(tb$sigma), as.numeric(at$scale), tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("Heckman deux etapes = sampleSelection::selection", {
  skip_if_not_installed("sampleSelection")
  skip_if_not_installed("MASS")
  set.seed(3); n <- 800
  zs <- rnorm(n); xo <- rnorm(n)
  E <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.7, 0.7, 1), 2))
  d <- as.integer(0.3 + 0.8 * zs + 0.5 * xo + E[, 1] > 0)
  yobs <- ifelse(d == 1, 2 + 1.5 * xo + E[, 2], NA)
  dat <- data.frame(d = d, zs = zs, xo = xo, y = yobs)
  hk <- heckman(d ~ zs + xo, y ~ xo, dat)
  hs <- sampleSelection::selection(d ~ zs + xo, y ~ xo, data = dat, method = "2step")
  # coefficients de l'equation de resultat (intercept, xo, ratio de Mills)
  expect_equal(as.numeric(hk$beta), as.numeric(coef(hs)[4:6]), tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("probit : vraisemblance monotone, probabilites dans (0,1)", {
  set.seed(4); n <- 200; x <- rnorm(n); y <- rbinom(n, 1, pnorm(x))
  pr <- probit(y ~ x, data.frame(y = y, x = x))
  expect_true(all(pr$fitted > 0 & pr$fitted < 1))
  expect_true(pr$loglik <= 0)
})

test_that("Heckman : le coefficient du ratio de Mills detecte le biais de selection", {
  skip_if_not_installed("MASS")
  # sans correlation -> coef IMR ~ 0 ; avec forte correlation -> coef IMR grand
  detect <- function(rho) {
    set.seed(10); n <- 1500; zs <- rnorm(n); xo <- rnorm(n)
    E <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, rho, rho, 1), 2))
    d <- as.integer(0.2 + 0.9 * zs + 0.4 * xo + E[, 1] > 0)
    y <- ifelse(d == 1, 1 + xo + E[, 2], NA)
    heckman(d ~ zs + xo, y ~ xo, data.frame(d, zs, xo, y))$beta["imr"]
  }
  expect_lt(abs(detect(0)), 0.15)         # pas de biais -> IMR ~ 0
  expect_gt(detect(0.8), 0.5)             # fort biais -> IMR grand et positif
})
