# Tests — Module 27 (processus gaussiens). Référence : DiceKriging (covariance gauss).

make_gp <- function(n = 15, seed = 1) {
  set.seed(seed)
  X <- matrix(sort(runif(n, 0, 5)), n, 1); y <- sin(X[, 1]) + 0.1 * rnorm(n)
  Xs <- matrix(seq(0, 5, length.out = 30), 30, 1)
  list(X = X, y = y, Xs = Xs)
}

test_that("moyenne et variance du GP = DiceKriging::km à 1e-8", {
  skip_if_not_installed("DiceKriging")
  d <- make_gp(); l <- 1; sf <- 1; sn <- 0.1
  pr <- gp_predict(gp_fit(d$X, d$y, l, sf, sn), d$Xs)
  km <- DiceKriging::km(design = as.data.frame(d$X), response = d$y, covtype = "gauss",
                        coef.trend = 0, coef.cov = l, coef.var = sf^2, nugget = sn^2,
                        control = list(trace = FALSE))
  dk <- DiceKriging::predict(km, newdata = as.data.frame(d$Xs), type = "SK", checkNames = FALSE)
  expect_lt(max(abs(pr$mean - dk$mean)), 1e-8)
  expect_lt(max(abs(pr$sd_obs^2 - dk$sd^2)), 1e-8)          # variance d'observation
})

test_that("moyenne a posteriori GP = régression ridge à noyau (Prop. 27.2)", {
  d <- make_gp(); l <- 1.2; sf <- 1; sn <- 0.2
  pr <- gp_predict(gp_fit(d$X, d$y, l, sf, sn), d$Xs)
  kr <- kernel_ridge(d$X, d$y, l, sf, lambda = sn^2)
  expect_lt(max(abs(pr$mean - kr(d$Xs))), 1e-8)
})

test_that("vraisemblance marginale = calcul direct (éq. 27.3)", {
  d <- make_gp(); l <- 0.8; sf <- 1.3; sn <- 0.15; n <- nrow(d$X)
  fit <- gp_fit(d$X, d$y, l, sf, sn)
  K <- rbf_kernel(d$X, d$X, l, sf^2) + sn^2 * diag(n)
  ll <- -0.5 * as.numeric(t(d$y) %*% solve(K, d$y)) -
    0.5 * as.numeric(determinant(K, logarithm = TRUE)$modulus) - 0.5 * n * log(2 * pi)
  expect_equal(fit$loglik, ll, tolerance = 1e-8)
})

test_that("interpolation : sans bruit, la moyenne passe par les points", {
  d <- make_gp()
  pr <- gp_predict(gp_fit(d$X, d$y, lengthscale = 0.3, sigma_f = 1, sigma_n = 1e-4), d$X)
  expect_lt(max(abs(pr$mean - d$y)), 1e-2)                  # interpole (bruit négligeable)
  expect_lt(max(pr$sd), 0.05)                               # incertitude ~ 0 aux données
})

test_that("optimisation de la vraisemblance marginale : loglik ne diminue pas", {
  d <- make_gp()
  ll0 <- gp_fit(d$X, d$y, 1, 1, 0.1)$loglik
  o <- gp_optimize(d$X, d$y, init = c(1, 1, 0.1))
  expect_gte(o$loglik, ll0 - 1e-6)                          # optim >= init
  expect_true(all(c(o$lengthscale, o$sigma_f, o$sigma_n) > 0))
})
