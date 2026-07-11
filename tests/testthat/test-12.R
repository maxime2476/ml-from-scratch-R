# Tests — Module 12 (MLP). From scratch : pas de référence externe.
# Le test central est la VÉRIFICATION DU GRADIENT par différences finies.

relerr <- function(a, b) {
  fa <- unlist(a); fb <- unlist(b)
  sqrt(sum((fa - fb)^2)) / (sqrt(sum(fa^2)) + sqrt(sum(fb^2)) + 1e-30)
}

test_that("VÉRIF. GRADIENT : backprop = différences finies (MSE, toutes activations)", {
  set.seed(1); n <- 25; d0 <- 3; d1 <- 6
  X <- matrix(rnorm(n * d0), n, d0); Y <- matrix(rnorm(n), n, 1)
  params <- list(W1 = matrix(rnorm(d0 * d1), d0, d1), b1 = rnorm(d1),
                 W2 = matrix(rnorm(d1), d1, 1), b2 = rnorm(1))
  for (act in c("tanh", "relu", "sigmoid")) {
    ga <- mlp_backward(params, X, Y, act, "mse")
    gn <- mlp_numgrad(params, X, Y, act, "mse")
    expect_lt(relerr(ga, gn), 1e-6)
  }
})

test_that("VÉRIF. GRADIENT : backprop = différences finies (log-loss)", {
  set.seed(2); n <- 25; d0 <- 4; d1 <- 5
  X <- matrix(rnorm(n * d0), n, d0); Y <- matrix(rbinom(n, 1, 0.5), n, 1)
  params <- list(W1 = matrix(rnorm(d0 * d1), d0, d1), b1 = rnorm(d1),
                 W2 = matrix(rnorm(d1), d1, 1), b2 = rnorm(1))
  for (act in c("tanh", "sigmoid")) {
    ga <- mlp_backward(params, X, Y, act, "logloss")
    gn <- mlp_numgrad(params, X, Y, act, "logloss")
    expect_lt(relerr(ga, gn), 1e-6)
  }
})

test_that("passe avant : dimensions correctes", {
  params <- list(W1 = matrix(0, 3, 5), b1 = rep(0, 5),
                 W2 = matrix(0, 5, 1), b2 = 0)
  fw <- mlp_forward(matrix(rnorm(10 * 3), 10, 3), params = params, activation = "tanh")
  expect_equal(dim(fw$Z1), c(10, 5))
  expect_equal(dim(fw$A1), c(10, 5))
  expect_equal(dim(fw$Z2), c(10, 1))
})

test_that("le MLP apprend le XOR (non linéairement séparable)", {
  set.seed(3)
  g <- expand.grid(x1 = c(0, 1), x2 = c(0, 1))
  X <- g[rep(1:4, each = 60), ] + matrix(rnorm(240 * 2, 0, 0.05), ncol = 2)
  y <- as.integer(xor(g$x1 == 1, g$x2 == 1))[rep(1:4, each = 60)]
  m <- mlp_fit(as.matrix(X), matrix(y, ncol = 1), hidden = 8, activation = "tanh",
               loss = "logloss", epochs = 400, lr = 0.1, batch = 32, seed = 5)
  expect_gt(mean(predict_mlp(m, as.matrix(X), "class") == y), 0.95)
})

test_that("le MLP apprend une régression non linéaire mieux qu'un modèle linéaire", {
  set.seed(4)
  x <- matrix(runif(300, -3, 3), 300, 1); y <- sin(2 * x) + rnorm(300, sd = 0.1)
  m <- mlp_fit(x, y, hidden = 15, activation = "tanh", loss = "mse",
               epochs = 500, lr = 0.03, batch = 32, seed = 6)
  mse_mlp <- mean((predict_mlp(m, x) - y)^2)
  mse_lin <- mean(lm(y ~ x)$residuals^2)
  expect_lt(mse_mlp, 0.5 * mse_lin)                    # bat nettement le linéaire
})

test_that("la perte d'entraînement décroît globalement", {
  set.seed(7)
  x <- matrix(runif(200, -2, 2), 200, 1); y <- x^2 + rnorm(200, sd = 0.1)
  m <- mlp_fit(x, y, hidden = 10, loss = "mse", epochs = 200, lr = 0.03, seed = 8)
  expect_lt(m$loss_hist[200], m$loss_hist[1])
  expect_lt(mean(tail(m$loss_hist, 20)), mean(head(m$loss_hist, 20)))
})

test_that("predict_mlp : types response/class cohérents (log-loss)", {
  set.seed(9)
  X <- matrix(rnorm(80 * 2), 80, 2); y <- matrix(rbinom(80, 1, 0.5), 80, 1)
  m <- mlp_fit(X, y, hidden = 5, loss = "logloss", epochs = 50, seed = 1)
  p <- predict_mlp(m, X, "response"); cl <- predict_mlp(m, X, "class")
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(cl %in% c(0, 1)))
})
