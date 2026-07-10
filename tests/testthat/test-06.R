# Tests de conformité — Module 6 (validation). Références : LOOCV force brute,
# hatvalues, AIC()/BIC(). Tolérance 1e-8.

make_val <- function(n = 70, seed = 8) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- 1 - 2 * d$x1 + 0.5 * d$x2 + rnorm(n)
  d
}

test_that("loocv_linear = LOOCV par force brute (n ajustements)", {
  d <- make_val()
  X <- model.matrix(y ~ x1 + x2 + x3, d)
  lc <- loocv_linear(X, d$y)
  n <- nrow(X)
  brute <- numeric(n)
  for (i in seq_len(n)) {
    b <- solve_ls_qr(X[-i, ], d$y[-i])$coefficients
    brute[i] <- (d$y[i] - sum(X[i, ] * b))^2
  }
  expect_equal(lc$cv, mean(brute), tolerance = 1e-8)
})

test_that("leviers de loocv_linear = hatvalues(lm)", {
  d <- make_val()
  X <- model.matrix(y ~ x1 + x2 + x3, d)
  lc <- loocv_linear(X, d$y)
  expect_equal(unname(lc$h), unname(hatvalues(lm(y ~ x1 + x2 + x3, d))), tolerance = 1e-8)
})

test_that("gcv_linear cohérent (entre erreur d'apprentissage et LOOCV)", {
  d <- make_val()
  X <- model.matrix(y ~ x1 + x2 + x3, d)
  fit <- solve_ls_qr(X, d$y)
  train_err <- mean(fit$residuals^2)
  gc <- gcv_linear(X, d$y)$gcv
  lc <- loocv_linear(X, d$y)$cv
  expect_gt(gc, train_err)          # GCV pénalise l'optimisme
  expect_lt(abs(gc - lc), 0.5 * lc)  # proche du LOOCV
})

test_that("info_criteria reproduit AIC()/BIC() pour lm", {
  d <- make_val()
  fit <- ols_fit(y ~ x1 + x2 + x3, d)
  ic <- info_criteria(fit)
  ref <- lm(y ~ x1 + x2 + x3, d)
  expect_equal(ic$aic, AIC(ref), tolerance = 1e-8)
  expect_equal(ic$bic, BIC(ref), tolerance = 1e-8)
})

test_that("info_criteria reproduit AIC()/BIC() pour glm (binomial et Poisson)", {
  d <- make_val()
  d$yb <- rbinom(nrow(d), 1, plogis(d$x1))
  d$yp <- rpois(nrow(d), exp(0.3 + 0.4 * d$x1))
  gb <- glm_irls(yb ~ x1 + x2, d, "binomial")
  gp <- glm_irls(yp ~ x1 + x2, d, "poisson")
  expect_equal(info_criteria(gb)$aic, AIC(glm(yb ~ x1 + x2, d, family = binomial)), tolerance = 1e-8)
  expect_equal(info_criteria(gb)$bic, BIC(glm(yb ~ x1 + x2, d, family = binomial)), tolerance = 1e-8)
  expect_equal(info_criteria(gp)$aic, AIC(glm(yp ~ x1 + x2, d, family = poisson)), tolerance = 1e-8)
  expect_equal(info_criteria(gp)$bic, BIC(glm(yp ~ x1 + x2, d, family = poisson)), tolerance = 1e-8)
})

test_that("kfold_cv est reproductible et proche du LOOCV", {
  d <- make_val()
  X <- model.matrix(y ~ x1 + x2 + x3, d)
  k1 <- kfold_cv(X, d$y, K = 10, seed = 42)
  k2 <- kfold_cv(X, d$y, K = 10, seed = 42)
  expect_equal(k1$cv, k2$cv)                       # même graine -> même résultat
  expect_lt(abs(k1$cv - loocv_linear(X, d$y)$cv), 0.5)
})

test_that("bias_variance_mc : irréductible + biais² + variance = EQM de prédiction", {
  set.seed(3)
  n <- 40; sigma2 <- 1
  Xtrain_x <- rnorm(n)
  x0 <- 0.7; f <- function(x) 1 + 2 * x; f0 <- f(x0)
  # OLS y ~ x ; prédiction en x0
  gen_y <- function() f(Xtrain_x) + rnorm(n, sd = sqrt(sigma2))
  fit_pred <- function(y) {
    b <- solve_ls_qr(cbind(1, Xtrain_x), y)$coefficients
    b[1] + b[2] * x0
  }
  bv <- bias_variance_mc(gen_y, fit_pred, f0, sigma2, R = 20000)
  # EQM de prédiction directe : E[(y0 - pred)^2], y0 = f0 + bruit indépendant
  set.seed(5); R <- 20000
  errs <- numeric(R)
  for (r in seq_len(R)) {
    pred <- fit_pred(gen_y())
    y0 <- f0 + rnorm(1, sd = sqrt(sigma2))
    errs[r] <- (y0 - pred)^2
  }
  expect_equal(bv$total, mean(errs), tolerance = 0.05 * bv$total)
})
