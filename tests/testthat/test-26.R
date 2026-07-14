# Tests — Module 26 (double descente). Phénomène vérifié ; interpolant = OLS (M0/lm).

make_dd <- function(n = 80, p = 5, seed = 1) {
  set.seed(seed)
  f0 <- function(X) sin(2 * X[, 1]) + 0.5 * X[, 2]^2 - X[, 3]
  Xtr <- matrix(rnorm(n * p), n, p); ytr <- f0(Xtr) + 0.3 * rnorm(n)
  Xte <- matrix(rnorm(1000 * p), 1000, p); yte <- f0(Xte) + 0.3 * rnorm(1000)
  list(Xtr = Xtr, ytr = ytr, Xte = Xte, yte = yte, n = n)
}

test_that("caractéristiques fixes entre train et test (mêmes poids)", {
  set.seed(1); X1 <- matrix(rnorm(20), 4, 5); X2 <- matrix(rnorm(20), 4, 5)
  # même seed -> même carte ; deux appels sur le MÊME X donnent le même résultat
  expect_equal(random_features(X1, 10, seed = 3), random_features(X1, 10, seed = 3))
  expect_false(isTRUE(all.equal(random_features(X1, 10, seed = 3),
                                random_features(X2, 10, seed = 3))))
})

test_that("interpolant de norme minimale = OLS quand D < n", {
  d <- make_dd(); Phi <- random_features(d$Xtr, 30, seed = 7)   # D=30 < n=80
  b_svd <- solve_ls_svd(Phi, d$ytr)$coefficients
  b_lm  <- unname(coef(lm(d$ytr ~ Phi - 1)))
  expect_lt(max(abs(b_svd - b_lm)), 1e-8)
})

test_that("D > n : interpolation (erreur d'apprentissage ~ 0)", {
  d <- make_dd()
  pred <- fit_rff(d$Xtr, d$ytr, D = 200, seed = 7, lambda = 0)
  expect_lt(mean((d$ytr - pred(d$Xtr))^2), 1e-10)
})

test_that("double descente : pic au seuil D~n puis seconde descente", {
  d <- make_dd()
  Ds <- c(10, 40, 78, 80, 82, 120, 400)
  cur <- double_descent_curve(d$Xtr, d$ytr, d$Xte, d$yte, Ds, gamma = 0.5, seed = 7)
  peak <- cur$test_mse[which.min(abs(cur$D - d$n))]
  expect_equal(which.max(cur$test_mse), which.min(abs(cur$D - d$n)))  # max au seuil
  expect_lt(cur$test_mse[cur$D == 400], peak)                        # redescend sous le pic
  expect_lt(cur$test_mse[cur$D == 400], cur$test_mse[cur$D == 40])   # bat le régime sous-param.
})

test_that("le ridge supprime le pic (régularisation explicite)", {
  d <- make_dd()
  Ds <- c(10, 40, 80, 120, 400)
  cur <- double_descent_curve(d$Xtr, d$ytr, d$Xte, d$yte, Ds, gamma = 0.5, seed = 7, lambda = 1)
  expect_false(which.max(cur$test_mse) == which.min(abs(cur$D - d$n)))  # plus de pic au seuil
})
