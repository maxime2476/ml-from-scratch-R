# Tests — Module 8 (CART). Référence qualitative : rpart (premiers splits,
# performance). L'égalité parfaite n'est pas garantie (cf. dérivation).

skip_if_not_installed("rpart")

make_tree_data <- function(n = 400, seed = 5) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$yc <- factor(ifelse(d$x1 + 0.5 * d$x2 + rnorm(n) > 0, "A", "B"))
  d$yr <- d$x1^2 + 0.5 * d$x2 + rnorm(n, sd = 0.3)
  d
}

test_that("mesures d'impureté : valeurs connues", {
  y <- factor(c("a", "a", "b", "b"))                  # 50/50
  expect_equal(impurity_gini(y), 0.5)                 # 1 - (0.25+0.25)
  expect_equal(impurity_entropy(y), log(2))           # -2*0.5*log(0.5)
  expect_equal(impurity_gini(factor(c("a", "a", "a"))), 0)   # pur
  expect_equal(impurity_variance(c(1, 3, 5)), mean((c(1,3,5) - 3)^2))
})

test_that("best_split trouve le seuil séparateur évident", {
  X <- matrix(c(1, 2, 3, 10, 11, 12), 6, 1)
  y <- factor(c("a", "a", "a", "b", "b", "b"))
  bs <- best_split(X, y, "class", "gini", min_leaf = 1, classes = c("a", "b"))
  expect_equal(bs$var, 1)
  expect_true(bs$val > 3 && bs$val < 10)              # sépare a de b
  expect_equal(bs$gain, 0.5)                          # Gini 0.5 -> 0
})

test_that("CART classification : premier split identique à rpart", {
  d <- make_tree_data()
  fit <- cart_fit(yc ~ x1 + x2 + x3, d, "class", min_split = 20, min_leaf = 7)
  rp <- rpart::rpart(yc ~ x1 + x2 + x3, d, method = "class",
                     control = rpart::rpart.control(minsplit = 20, minbucket = 7, cp = 0, maxdepth = 30))
  expect_equal(fit$tree$var, as.character(rp$frame$var[1]))          # même variable
  expect_equal(fit$tree$val, unname(rp$splits[1, "index"]), tolerance = 1e-6)  # même seuil
})

test_that("CART régression : premier split et performance ~ rpart", {
  d <- make_tree_data()
  fit <- cart_fit(yr ~ x1 + x2 + x3, d, "anova", min_split = 20, min_leaf = 10)
  rp <- rpart::rpart(yr ~ x1 + x2 + x3, d, method = "anova",
                     control = rpart::rpart.control(minsplit = 20, minbucket = 10, cp = 0, maxdepth = 30))
  expect_equal(fit$tree$var, as.character(rp$frame$var[1]))
  expect_equal(fit$tree$val, unname(rp$splits[1, "index"]), tolerance = 1e-6)
  mse_mine <- mean((predict_cart(fit, d) - d$yr)^2)
  mse_rp   <- mean((predict(rp, d) - d$yr)^2)
  expect_lt(abs(mse_mine - mse_rp), 0.05 * mse_rp)
})

test_that("prédiction cohérente : classes valides, accuracy raisonnable", {
  d <- make_tree_data()
  fit <- cart_fit(yc ~ x1 + x2 + x3, d, "class")
  pred <- predict_cart(fit, d)
  expect_true(all(levels(pred) == c("A", "B")))
  expect_gt(mean(pred == d$yc), 0.75)
})

test_that("cost_complexity_prune : alpha=0 garde tout, alpha grand -> racine", {
  d <- make_tree_data()
  fit <- cart_fit(yr ~ x1 + x2 + x3, d, "anova", min_split = 10, min_leaf = 5)
  full <- n_leaves(fit)
  expect_equal(n_leaves(cost_complexity_prune(fit, 0)), full)        # rien élagué
  expect_lt(n_leaves(cost_complexity_prune(fit, 0.05)), full)        # élagage
  expect_equal(n_leaves(cost_complexity_prune(fit, 1e6)), 1L)        # réduit à la racine
})

test_that("contre-exemple XOR : greedy myope en profondeur 1, résolu en profondeur 2", {
  set.seed(1)
  m <- 100
  g <- expand.grid(x1 = c(0, 1), x2 = c(0, 1))
  d <- g[rep(seq_len(4), each = m), ]
  d$y <- factor(as.integer(xor(d$x1 == 1, d$x2 == 1)))
  # profondeur 1 : aucun split à gain > 0 -> feuille -> ~50 %
  f1 <- cart_fit(y ~ x1 + x2, d, "class", max_depth = 1, min_split = 2, min_leaf = 1, min_gain = -1)
  # profondeur 2 en acceptant le split à gain nul (min_gain < 0) -> parfait
  f2 <- cart_fit(y ~ x1 + x2, d, "class", max_depth = 2, min_split = 2, min_leaf = 1, min_gain = -1)
  expect_equal(mean(predict_cart(f1, d) == d$y), 0.5, tolerance = 1e-8)
  expect_equal(mean(predict_cart(f2, d) == d$y), 1.0, tolerance = 1e-8)
})
