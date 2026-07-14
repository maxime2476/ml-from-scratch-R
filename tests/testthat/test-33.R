# Tests — Module 33 (MCMC). References : moments cibles analytiques, OLS, coda.

test_that("Metropolis-Hastings echantillonne la cible (moments recuperes)", {
  set.seed(1)
  r <- metropolis_hastings(function(x) dnorm(x, 3, 2, log = TRUE), 0, 2.5, 30000)
  ch <- r$chain[5001:30000, 1]
  expect_lt(abs(mean(ch) - 3), 0.1)          # moyenne cible = 3
  expect_lt(abs(sd(ch) - 2), 0.1)            # ecart-type cible = 2
  expect_true(r$accept_rate > 0.2 && r$accept_rate < 0.9)
})

test_that("MH sur une cible Beta (bornee) : moments corrects", {
  set.seed(2)
  a <- 2; b <- 5
  lt <- function(x) if (x > 0 && x < 1) (a - 1) * log(x) + (b - 1) * log(1 - x) else -Inf
  ch <- metropolis_hastings(lt, 0.3, 0.1, 40000)$chain[10001:40000, 1]
  expect_lt(abs(mean(ch) - a / (a + b)), 0.02)              # E = a/(a+b)
})

test_that("Gibbs (regression bayesienne) : moyenne a posteriori = OLS", {
  set.seed(3); n <- 120; X <- cbind(1, rnorm(n), rnorm(n))
  y <- as.numeric(X %*% c(1, 2, -1)) + rnorm(n)
  g <- gibbs_linreg(X, y, n_iter = 6000, burn = 1000)
  ols <- as.numeric(solve(crossprod(X), crossprod(X, y)))
  expect_lt(max(abs(colMeans(g$beta) - ols)), 0.05)
  expect_lt(abs(mean(g$sigma2) - sum((y - X %*% ols)^2) / (n - 3)) / mean(g$sigma2), 0.15)
})

test_that("Gelman-Rubin ~ 1 a convergence, > 1 si chaines disjointes", {
  set.seed(4)
  lt <- function(x) dnorm(x, 0, 1, log = TRUE)
  c1 <- metropolis_hastings(lt, 0, 1.5, 8000)$chain[2001:8000, 1]
  c2 <- metropolis_hastings(lt, 0, 1.5, 8000)$chain[2001:8000, 1]
  expect_lt(gelman_rubin(list(c1, c2)), 1.05)               # converge
  # chaines bloquees sur des modes differents -> R-hat grand
  stuck <- gelman_rubin(list(rnorm(2000, -5), rnorm(2000, 5)))
  expect_gt(stuck, 1.5)
})

test_that("Gelman-Rubin et ESS ~ coda", {
  skip_if_not_installed("coda")
  set.seed(5); lt <- function(x) dnorm(x, 3, 2, log = TRUE)
  c1 <- metropolis_hastings(lt, 0, 2.5, 10000)$chain[2001:10000, 1]
  c2 <- metropolis_hastings(lt, 6, 2.5, 10000)$chain[2001:10000, 1]
  rh <- gelman_rubin(list(c1, c2))
  rc <- coda::gelman.diag(coda::mcmc.list(coda::mcmc(c1), coda::mcmc(c2)))$psrf[1]
  expect_lt(abs(rh - rc), 0.02)
  expect_lt(abs(ess(c1) - coda::effectiveSize(c1)) / coda::effectiveSize(c1), 0.25)
})
