# Tests — Module 9 (bagging / forêts aléatoires). Référence qualitative :
# randomForest (aléa -> pas d'égalité exacte).

skip_if_not_installed("randomForest")

make_rf_data <- function(n = 400, seed = 7) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n))
  d$yc <- factor(ifelse(d$x1 + d$x2 - 0.5 * d$x3 + rnorm(n) > 0, "A", "B"))
  d$yr <- d$x1^2 + d$x2 - 0.5 * d$x3 * d$x4 + rnorm(n, sd = 0.4)
  d
}

test_that("fraction OOB moyenne ~ 1/e (Prop. 9.2)", {
  fit <- bagging_fit(yc ~ x1 + x2 + x3 + x4, make_rf_data(), "class", B = 200, seed = 1)
  oob_frac <- mean(vapply(fit$oob_idx, length, integer(1))) / fit$n
  expect_equal(oob_frac, 1 / exp(1), tolerance = 0.02)
})

test_that("forêt de classification ~ randomForest (OOB et accuracy test)", {
  d <- make_rf_data(); set.seed(1)
  te <- make_rf_data(500, seed = 99)
  fit <- random_forest_fit(yc ~ x1 + x2 + x3 + x4, d, "class", B = 200, seed = 1)
  rf <- randomForest::randomForest(yc ~ x1 + x2 + x3 + x4, d, ntree = 200)
  expect_lt(abs(fit$oob_error - rf$err.rate[200, "OOB"]), 0.05)   # OOB proche
  acc_mine <- mean(predict_forest(fit, te) == te$yc)
  acc_rf   <- mean(predict(rf, te) == te$yc)
  expect_lt(abs(acc_mine - acc_rf), 0.05)
})

test_that("forêt de régression ~ randomForest (OOB MSE)", {
  d <- make_rf_data()
  fr <- random_forest_fit(yr ~ x1 + x2 + x3 + x4, d, "anova", B = 200, seed = 2)
  rfr <- randomForest::randomForest(yr ~ x1 + x2 + x3 + x4, d, ntree = 200)
  expect_lt(abs(fr$oob_error - rfr$mse[200]) / rfr$mse[200], 0.15)  # ~15 %
})

test_that("le bagging réduit la variance vs un arbre seul", {
  # Variance de la prédiction en des points de test, sur plusieurs jeux d'apprentissage.
  set.seed(3)
  x0 <- data.frame(x1 = 0.5, x2 = -0.5, x3 = 0.2, x4 = 0.1)
  R <- 40
  p_tree <- numeric(R); p_bag <- numeric(R)
  for (r in seq_len(R)) {
    d <- make_rf_data(200, seed = r)
    tr <- cart_fit(yr ~ x1 + x2 + x3 + x4, d, "anova", max_depth = 20, min_split = 5, min_leaf = 1)
    bg <- bagging_fit(yr ~ x1 + x2 + x3 + x4, d, "anova", B = 50, mtry = 4, seed = r)
    p_tree[r] <- predict_cart(tr, x0)
    p_bag[r]  <- predict_forest(bg, x0)
  }
  expect_lt(var(p_bag), var(p_tree))            # variance réduite par agrégation
})

test_that("predict_forest : sorties valides (classes / numériques)", {
  d <- make_rf_data()
  fc <- random_forest_fit(yc ~ x1 + x2 + x3 + x4, d, "class", B = 50, seed = 1)
  fr <- random_forest_fit(yr ~ x1 + x2 + x3 + x4, d, "anova", B = 50, seed = 1)
  pc <- predict_forest(fc, d[1:10, ]); pr <- predict_forest(fr, d[1:10, ])
  expect_true(all(levels(pc) == c("A", "B")))
  expect_true(is.numeric(pr) && length(pr) == 10)
})
