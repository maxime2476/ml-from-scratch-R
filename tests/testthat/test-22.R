# Tests — Module 22 (lasso débiaisé). Références : hdm, couverture MC.

make_hd <- function(n = 100, p = 200, s = 5, seed = 1) {
  set.seed(seed)
  beta <- c(rep(1.5, s), rep(0, p - s))
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% beta + rnorm(n))
  list(X = X, y = y, beta = beta)
}

test_that("debiased_lasso : IC couvrent les vraies valeurs (actif et nul)", {
  d <- make_hd()
  db <- debiased_lasso(d$X, d$y, targets = c(1, 6))   # b1=1.5, b6=0
  expect_true(db$lower[1] <= 1.5 && 1.5 <= db$upper[1])
  expect_true(db$lower[2] <= 0 && 0 <= db$upper[2])
})

test_that("estimations ponctuelles ~ hdm::rlassoEffects", {
  skip_if_not_installed("hdm")
  d <- make_hd(seed = 2)
  db <- debiased_lasso(d$X, d$y, targets = 1:6)
  he <- hdm::rlassoEffects(d$X, d$y, index = 1:6)
  # même ordre de grandeur, corrélés (les deux estiment les mêmes effets)
  expect_gt(cor(db$estimate, he$coefficients), 0.8)
  expect_lt(max(abs(db$estimate - he$coefficients)), 0.3)
})

test_that("couverture Monte Carlo ~ 0.95 (actif et nul) en p > n", {
  R <- 300; ca <- cn <- 0
  for (r in seq_len(R)) {
    d <- make_hd(seed = 100 + r)
    db <- debiased_lasso(d$X, d$y, targets = c(1, 6))
    ca <- ca + (db$lower[1] <= 1.5 && 1.5 <= db$upper[1])
    cn <- cn + (db$lower[2] <= 0 && 0 <= db$upper[2])
  }
  expect_gt(ca / R, 0.90)          # couverture valide (vs naïf qui échoue, M14)
  expect_gt(cn / R, 0.90)
})
