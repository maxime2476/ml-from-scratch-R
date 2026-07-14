# Tests — Module 34 (methode delta, tests multiples). References : car, p.adjust.

test_that("methode delta = car::deltaMethod (ratio de coefficients)", {
  skip_if_not_installed("car")
  set.seed(1); n <- 200; x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 1 + 2 * x1 + 3 * x2 + rnorm(n); fit <- lm(y ~ x1 + x2)
  dm <- delta_method(coef(fit), vcov(fit), function(b) b[2] / b[3])
  dc <- car::deltaMethod(fit, "x1/x2")
  expect_equal(dm$estimate, as.numeric(dc$Estimate), tolerance = 1e-5)
  expect_equal(dm$se, as.numeric(dc$SE), tolerance = 1e-4)
})

test_that("methode delta : cas scalaire connu (exp(theta))", {
  # Var(exp(theta)) = exp(theta)^2 * Var(theta) (delta exact au 1er ordre)
  dm <- delta_method(c(0.5), matrix(0.04), function(b) exp(b[1]))
  expect_equal(dm$se, exp(0.5) * sqrt(0.04), tolerance = 1e-5)
})

test_that("Bonferroni et BH = stats::p.adjust", {
  set.seed(2); p <- c(runif(15, 0, 0.01), runif(85))
  expect_equal(p_adjust_bonferroni(p), p.adjust(p, "bonferroni"), tolerance = 1e-12)
  expect_equal(p_adjust_bh(p), p.adjust(p, "BH"), tolerance = 1e-12)
})

test_that("BH est plus puissant que Bonferroni (plus de rejets)", {
  set.seed(3); p <- c(runif(30, 0, 0.005), runif(70))       # 30 vrais positifs
  n_bonf <- sum(p_adjust_bonferroni(p) < 0.05)
  n_bh <- sum(p_adjust_bh(p) < 0.05)
  expect_gte(n_bh, n_bonf)                                   # BH rejette au moins autant
})
