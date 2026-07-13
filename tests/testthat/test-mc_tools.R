# Tests — outils de rigueur Monte Carlo.

test_that("mc_se = sd/sqrt(R)", {
  set.seed(1); x <- rnorm(500)
  expect_equal(mc_se(x), sd(x) / sqrt(500))
})

test_that("mc_summary : biais, RMSE et leurs erreurs MC corrects", {
  set.seed(1); est <- rnorm(2000, mean = 2.1, sd = 0.3); truth <- 2
  s <- mc_summary(est, truth)
  expect_equal(s$bias, mean(est) - 2, tolerance = 1e-12)
  expect_equal(s$bias_se, sd(est) / sqrt(2000), tolerance = 1e-12)
  expect_equal(s$rmse, sqrt(mean((est - 2)^2)), tolerance = 1e-12)
  expect_gt(s$rmse_se, 0)
})

test_that("coverage_mc : proportion, erreur MC binomiale, test du nominal", {
  set.seed(2); cov <- rbinom(1000, 1, 0.95) == 1
  cm <- coverage_mc(cov, nominal = 0.95)
  expect_equal(cm$coverage, mean(cov))
  expect_equal(cm$se, sqrt(mean(cov) * (1 - mean(cov)) / 1000), tolerance = 1e-12)
  expect_true(cm$nominal_ok)                          # ~0.95 : IC contient 0.95
  # couverture nettement basse -> nominal rejeté
  bad <- coverage_mc(rbinom(2000, 1, 0.80) == 1, nominal = 0.95)
  expect_false(bad$nominal_ok)
})

test_that("reject_mc : taux et erreur MC", {
  set.seed(3); rej <- rbinom(4000, 1, 0.05) == 1
  rm <- reject_mc(rej, nominal = 0.05)
  expect_equal(rm$rate, mean(rej))
  expect_true(rm$nominal_ok)                          # taille ~ nominale
})

test_that("convergence_study : OLS est sqrt(n)-consistant (biais->0, sqrt(n)*sd stable)", {
  # sim_fn : coefficient OLS de x sur un DGP simple, theta = 2
  sim_ols <- function(n) {
    x <- rnorm(n); y <- 1 + 2 * x + rnorm(n)
    coef(lm(y ~ x))[2]
  }
  conv <- convergence_study(sim_ols, ns = c(50, 200, 800), R = 400, truth = 2, seed = 1)
  expect_true(all(abs(conv$bias) < 3 * conv$bias_se + 0.02))   # biais ~ 0 partout
  # sqrt(n)*sd ~ constante (ecart-type asymptotique) : faible variation relative
  expect_lt(sd(conv$sqrtn_sd) / mean(conv$sqrtn_sd), 0.15)
  # pente log-log du RMSE ~ -0.5
  expect_lt(abs(rmse_rate(conv) + 0.5), 0.12)
})
