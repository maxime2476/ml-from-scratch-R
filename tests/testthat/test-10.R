# Tests — Module 10 (gradient boosting). Math exacte + validation qualitative gbm.

make_boost_data <- function(n = 500, seed = 4) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$yr <- sin(2 * d$x1) + d$x2 - 0.5 * d$x3 + rnorm(n, sd = 0.4)
  d$yb <- rbinom(n, 1, plogis(1.2 * d$x1 - 0.8 * d$x2))
  d
}

test_that("L2 : une itération = F0 + arbre ajusté aux résidus (éq. 10.3)", {
  d <- make_boost_data(200)
  fit <- gradient_boost(yr ~ x1 + x2 + x3, d, "l2", M = 1, nu = 1,
                        max_depth = 3, min_leaf = 10)
  F0 <- mean(d$yr)
  tr <- cart_fit(r ~ x1 + x2 + x3, cbind(d, r = d$yr - F0), "anova",
                 max_depth = 3, min_leaf = 10)
  expect_equal(predict_boost(fit, d, "link"), F0 + predict_cart(tr, d), tolerance = 1e-10)
})

test_that("Newton : poids de feuille = -sum(g)/(sum(h)+lambda) (éq. 10.7, log-loss)", {
  d <- make_boost_data(300)
  # arbre-racine (profondeur 0) : une seule feuille sur tout l'échantillon
  fit <- gradient_boost(yb ~ x1 + x2 + x3, d, "logloss", M = 1, nu = 1,
                        max_depth = 0, newton = TRUE, lambda = 0)
  p0 <- mean(d$yb)                                   # p = sigma(F0) = ybar
  w_manuel <- sum(d$yb - p0) / sum(p0 * (1 - p0))     # = sum(y-p)/sum(p(1-p))
  F0 <- log(p0 / (1 - p0))
  expect_equal(unique(predict_boost(fit, d, "link")), F0 + w_manuel, tolerance = 1e-8)
})

test_that("log-loss : la perte d'entraînement décroît de façon monotone", {
  d <- make_boost_data()
  fit <- gradient_boost(yb ~ x1 + x2 + x3, d, "logloss", M = 100, nu = 0.1,
                        max_depth = 3, min_leaf = 10)
  lp <- boost_loss_path(fit, d, d$yb)
  expect_true(all(diff(lp) <= 1e-8))                 # décroissance
  expect_lt(lp[100], lp[1])
})

test_that("taux d'apprentissage : petit nu ajuste moins vite (à M fixe)", {
  d <- make_boost_data()
  f_fast <- gradient_boost(yr ~ x1 + x2 + x3, d, "l2", M = 50, nu = 0.3, max_depth = 3)
  f_slow <- gradient_boost(yr ~ x1 + x2 + x3, d, "l2", M = 50, nu = 0.02, max_depth = 3)
  l_fast <- boost_loss_path(f_fast, d, d$yr)[50]
  l_slow <- boost_loss_path(f_slow, d, d$yr)[50]
  expect_gt(l_slow, l_fast)                           # nu petit -> perte train plus haute
})

test_that("predict_boost : types response/link/class cohérents", {
  d <- make_boost_data(200)
  fb <- gradient_boost(yb ~ x1 + x2 + x3, d, "logloss", M = 50, nu = 0.1)
  p <- predict_boost(fb, d, "response"); cl <- predict_boost(fb, d, "class")
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(cl %in% c(0L, 1L)))
  expect_equal(cl, as.integer(p >= 0.5))
})

test_that("L2 : performance et trajectoire de perte ~ gbm", {
  skip_if_not_installed("gbm")
  d <- make_boost_data(); set.seed(5)
  te <- make_boost_data(1000, seed = 11)
  fit <- gradient_boost(yr ~ x1 + x2 + x3, d, "l2", M = 200, nu = 0.1,
                        max_depth = 3, min_leaf = 10)
  gb <- gbm::gbm(yr ~ x1 + x2 + x3, data = d, distribution = "gaussian", n.trees = 200,
                 shrinkage = 0.1, interaction.depth = 3, n.minobsinnode = 10, bag.fraction = 1)
  mse_mine <- mean((predict_boost(fit, te) - te$yr)^2)
  mse_gbm  <- mean((gbm::predict.gbm(gb, te, n.trees = 200) - te$yr)^2)
  expect_lt(abs(mse_mine - mse_gbm) / mse_gbm, 0.20)                 # ~20 %
  expect_gt(cor(boost_loss_path(fit, d, d$yr), gb$train.error), 0.95) # trajectoires
})
