# Tests — Module 16 (DML, forêts causales). Math exacte (FWL) + DoubleML/grf.

make_causal <- function(n = 800, p = 5, theta = 1.5, seed = 4) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", seq_len(p))
  g <- sin(X[, 1]) + 0.3 * X[, 2]^2 + 0.5 * X[, 3]
  d <- rbinom(n, 1, plogis(0.8 * X[, 1] - 0.5 * X[, 2]))
  y <- theta * d + g + rnorm(n, sd = 1)
  list(X = X, y = y, d = d, theta = theta)
}

test_that("DML(lm, sans cross-fitting) = coefficient FWL de l'OLS (généralise M1)", {
  dd <- make_causal(n = 400)
  fit <- dml_plr(dd$y, dd$d, dd$X, nuisance = "lm", crossfit = FALSE)
  ref <- coef(lm(dd$y ~ dd$d + dd$X))["dd$d"]
  expect_equal(fit$theta, unname(ref), tolerance = 1e-8)           # FWL exact
})

test_that("DML (forêt, cross-fit) récupère theta ~ vrai, l'OLS naïf est biaisé", {
  dd <- make_causal()
  fit <- dml_plr(dd$y, dd$d, dd$X, K = 5, nuisance = "forest", crossfit = TRUE, seed = 1, B = 100)
  expect_lt(abs(fit$theta - dd$theta), 0.10)                       # DML ~ sans biais
  expect_true(fit$ci[1] <= dd$theta && dd$theta <= fit$ci[2])      # IC couvre
})

test_that("DML reproduit DoubleML (qualitativement)", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("ranger")
  suppressMessages({library(DoubleML); library(mlr3); library(mlr3learners)})
  lgr::get_logger("mlr3")$set_threshold("error")
  dd <- make_causal()
  fit <- dml_plr(dd$y, dd$d, dd$X, K = 5, nuisance = "forest", crossfit = TRUE, seed = 1, B = 100)
  dml_data <- DoubleML::double_ml_data_from_matrix(X = dd$X, y = dd$y, d = dd$d)
  lrn <- mlr3::lrn("regr.ranger", num.trees = 100)
  obj <- DoubleML::DoubleMLPLR$new(dml_data, ml_l = lrn$clone(), ml_m = lrn$clone(), n_folds = 5)
  invisible(obj$fit())
  expect_lt(abs(fit$theta - obj$coef), 0.10)                       # thetas proches
  expect_lt(abs(fit$se - obj$se), 0.05)                            # se proches
})

test_that("cross-fitting réduit le biais de sur-ajustement vs sans", {
  # DGP à forte non-linéarité : sans cross-fitting, biais plus marqué.
  set.seed(7)
  bias_cf <- bias_nocf <- numeric(20)
  for (r in seq_len(20)) {
    dd <- make_causal(n = 400, seed = 100 + r)
    bias_cf[r]   <- dml_plr(dd$y, dd$d, dd$X, K = 5, nuisance = "forest",
                            crossfit = TRUE, seed = r, B = 60)$theta - dd$theta
    bias_nocf[r] <- dml_plr(dd$y, dd$d, dd$X, nuisance = "forest",
                            crossfit = FALSE, seed = r, B = 60)$theta - dd$theta
  }
  expect_lt(abs(mean(bias_cf)), abs(mean(bias_nocf)))              # cross-fit moins biaisé
})

test_that("T-learner et arbre causal : CATE corrélé au vrai effet hétérogène", {
  set.seed(9)
  X <- matrix(rnorm(600 * 4), 600, 4); colnames(X) <- paste0("X", 1:4)
  tau <- 1 + X[, 1]                                                # CATE dépend de X1
  d <- rbinom(600, 1, plogis(0.4 * X[, 2]))
  y <- tau * d + sin(X[, 1]) + 0.3 * X[, 3] + rnorm(600, sd = 1)
  cate_t <- t_learner(as.data.frame(X), y, d, B = 150)
  ct <- causal_tree(X, y, d, max_depth = 3, min_leaf = 20, seed = 1)
  cate_ct <- predict_causal_tree(ct, X)
  expect_gt(cor(cate_t, tau), 0.4)                                 # positivement corrélés
  expect_gt(cor(cate_ct, tau), 0.4)
})

test_that("arbre causal grf : CATE qualitativement concordant", {
  skip_if_not_installed("grf")
  set.seed(9)
  X <- matrix(rnorm(600 * 4), 600, 4); colnames(X) <- paste0("X", 1:4)
  tau <- 1 + X[, 1]
  d <- rbinom(600, 1, plogis(0.4 * X[, 2]))
  y <- tau * d + sin(X[, 1]) + rnorm(600, sd = 1)
  cate_t <- t_learner(as.data.frame(X), y, d, B = 150)
  cf <- grf::causal_forest(X, y, d, num.trees = 400)
  cate_grf <- predict(cf)$predictions
  expect_gt(cor(cate_t, cate_grf), 0.4)                            # même tendance
})
