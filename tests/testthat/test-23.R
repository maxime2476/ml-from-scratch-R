# Tests — Module 23 (analyse de sensibilité OVB). Référence : sensemakr.

make_sens <- function(n = 300, seed = 3) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$D <- 0.5 * d$x1 + rnorm(n)
  d$y <- 1 + 2 * d$D + 1.5 * d$x1 - d$x2 + rnorm(n)
  d
}

test_that("R² partiel = t²/(t²+df)", {
  expect_equal(partial_r2(4, 100), 16 / 116)
  expect_equal(partial_r2(0, 50), 0)
})

test_that("robustness value dans [0,1] et décroît avec df", {
  rv1 <- robustness_value(5, 50)
  rv2 <- robustness_value(5, 500)   # même t, plus de df -> effet plus fragile
  expect_true(rv1 >= 0 && rv1 <= 1)
  expect_gt(rv1, rv2)
})

test_that("estimation ajustée : biais nul si r2=0, réduit l'effet sinon", {
  a0 <- adjusted_estimate(2, 0.1, 100, 0, 0)
  expect_equal(a0$estimate, 2); expect_equal(a0$bias, 0)
  a1 <- adjusted_estimate(2, 0.1, 100, 0.2, 0.2)
  expect_lt(a1$estimate, 2)          # réduit vers 0
  expect_gt(a1$se, 0.1 * sqrt((1 - 0.2) / (1 - 0.2)) * 0)  # se ajustée définie
})

test_that("R² partiel, RV et estimation ajustée = sensemakr", {
  skip_if_not_installed("sensemakr")
  d <- make_sens()
  fit <- ols_fit(y ~ D + x1 + x2, d)
  lmf <- lm(y ~ D + x1 + x2, d)
  s <- sensitivity_ols(fit, "D")
  sm <- sensemakr::sensemakr(model = lmf, treatment = "D")$sensitivity_stats
  expect_equal(s$r2yd, as.numeric(sm$r2yd), tolerance = 1e-6)
  expect_equal(s$rv_q, as.numeric(sm$rv_q), tolerance = 1e-6)
  expect_equal(s$rv_qa, as.numeric(sm$rv_qa), tolerance = 1e-6)

  aj <- adjusted_estimate(s$estimate, s$se, s$df, 0.1, 0.1)
  expect_equal(aj$estimate,
               as.numeric(sensemakr::adjusted_estimate(estimate = s$estimate, se = s$se,
                                            dof = s$df, r2dz.x = 0.1, r2yz.dx = 0.1)),
               tolerance = 1e-6)
  expect_equal(aj$se,
               as.numeric(sensemakr::adjusted_se(se = s$se, dof = s$df,
                                      r2dz.x = 0.1, r2yz.dx = 0.1)),
               tolerance = 1e-6)
})
