# Tests — Module 15 (interprétabilité). Math exacte (forme fermée linéaire) +
# comparaison à iml.

make_interp_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n))
  X
}

test_that("Shapley exact = forme fermée linéaire beta_j (x_j - E[x_j]) (Prop. 15.3)", {
  X <- make_interp_data()
  beta <- c(2, -1.5, 0.8, 0)
  pf <- function(D) as.numeric(as.matrix(D[, c("x1", "x2", "x3", "x4")]) %*% beta)
  xstar <- X[1, ]
  phi <- shapley_exact(pf, xstar, X)
  closed <- beta * (as.numeric(xstar) - colMeans(X))
  expect_equal(unname(phi), unname(closed), tolerance = 1e-8)
})

test_that("Shapley : exactitude locale (efficience) sum(phi) = f(x) - E[f]", {
  X <- make_interp_data()
  pf <- function(D) 1 + as.matrix(D) %*% c(2, -1, 0.5, 0.3) - 0.4 * D$x1 * D$x2
  phi <- shapley_exact(function(D) as.numeric(pf(D)), X[3, ], X)
  expect_equal(sum(phi), as.numeric(pf(X[3, ])) - mean(as.numeric(pf(X))), tolerance = 1e-8)
})

test_that("Shapley par permutation converge vers l'exact", {
  X <- make_interp_data(200)
  beta <- c(1.5, -1, 0.7, 0.2)
  pf <- function(D) as.numeric(as.matrix(D[, 1:4]) %*% beta)
  phi_ex <- shapley_exact(pf, X[1, ], X)
  phi_pm <- shapley_permutation(pf, X[1, ], X, n_samples = 8000, seed = 1)
  expect_lt(max(abs(phi_pm - phi_ex)), 0.03)
})

test_that("PDP d'un modèle linéaire a pour pente beta_j", {
  X <- make_interp_data()
  pf <- function(D) as.numeric(as.matrix(D[, 1:4]) %*% c(2, -1.5, 0.8, 0))
  pd <- pdp(pf, X, "x1")
  expect_equal(unname(coef(lm(pdp ~ grid, pd))[2]), 2, tolerance = 1e-8)
})

test_that("ICE : le PDP est la moyenne des courbes ICE", {
  X <- make_interp_data(100)
  pf <- function(D) as.numeric(as.matrix(D[, 1:4]) %*% c(1, 2, -1, 0.5) + 0.3 * D$x1 * D$x3)
  ic <- ice(pf, X, "x1", grid_size = 15)
  expect_equal(ic$pdp, colMeans(ic$ice), tolerance = 1e-12)
  # comparaison au PDP direct
  pd <- pdp(pf, X, "x1", grid = ic$grid)
  expect_equal(ic$pdp, pd$pdp, tolerance = 1e-10)
})

test_that("importance par permutation : variable nulle ~ 0, ordre selon |beta|", {
  X <- make_interp_data()
  beta <- c(2, -1.5, 0.8, 0)
  pf <- function(D) as.numeric(as.matrix(D[, 1:4]) %*% beta)
  y <- pf(X) + rnorm(nrow(X), sd = 0.5)
  imp <- permutation_importance(pf, X, y, seed = 1)
  expect_lt(abs(imp["x4"]), 1e-8)                    # variable inutilisée
  expect_gt(imp["x1"], imp["x2"])                    # |2| > |1.5|
  expect_gt(imp["x2"], imp["x3"])                    # |1.5| > |0.8|
})

test_that("PDP et importance permutation reproduisent iml", {
  skip_if_not_installed("iml")
  X <- make_interp_data(200)
  y <- 2 * X$x1 - 1.5 * X$x2 + 0.5 * X$x1 * X$x2 + rnorm(200, sd = 0.3)
  mod <- lm(y ~ x1 * x2 + x3 + x4, data = cbind(X, y = y))
  pf <- function(D) as.numeric(predict(mod, newdata = as.data.frame(D)))
  pred_obj <- iml::Predictor$new(model = mod, data = X, y = y,
    predict.function = function(m, newdata) predict(m, newdata = as.data.frame(newdata)))
  grid <- seq(min(X$x1), max(X$x1), length.out = 20)
  mine <- pdp(pf, X, "x1", grid = grid)
  fe <- iml::FeatureEffect$new(pred_obj, feature = "x1", method = "pdp", grid.points = grid)
  expect_equal(mine$pdp, fe$results$.value, tolerance = 1e-8)           # PDP exact
  imp <- permutation_importance(pf, X, y, loss = function(yh, y) sqrt(mean((yh - y)^2)),
                                n_repeat = 30, seed = 1)
  fi <- iml::FeatureImp$new(pred_obj, loss = "rmse", n.repetitions = 30, compare = "difference")
  iml_imp <- fi$results$importance[match(names(imp), fi$results$feature)]
  expect_gt(cor(imp, iml_imp), 0.95)                                    # importances corrélées
})
