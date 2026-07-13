# Tests — Module 19 (prédiction conforme). La garantie (Th. 19.1) est
# vérifiée par simulation ; pas de référence externe.

# Modèle de base : OLS.
ols_fitfn  <- function(X, y) solve_ls_qr(cbind(1, X), y)$coefficients
ols_predfn <- function(b, X) as.numeric(cbind(1, X) %*% b)

test_that("conformal_quantile : bon ordre statistique et +Inf si rang > n", {
  s <- c(3, 1, 2, 5, 4)
  # n=5, alpha=0.1 : k = ceiling(0.9*6) = 6 > 5 -> Inf
  expect_equal(conformal_quantile(s, 0.1), Inf)
  # alpha=0.3 : k = ceiling(0.7*6) = 5 -> 5e plus petit = 5
  expect_equal(conformal_quantile(s, 0.3), 5)
  # alpha=0.5 : k = ceiling(0.5*6) = 3 -> 3e = 3
  expect_equal(conformal_quantile(s, 0.5), 3)
})

coverage_over_reps <- function(gen_xy, alpha = 0.1, n = 200, R = 1500, seed = 1) {
  set.seed(seed); cov <- logical(R)
  for (r in seq_len(R)) {
    d <- gen_xy(n + 1)
    tr <- 1:(n / 2); cal <- (n / 2 + 1):n; te <- n + 1
    cf <- conformal_split(d$X[tr, , drop = FALSE], d$y[tr],
                          d$X[cal, , drop = FALSE], d$y[cal],
                          d$X[te, , drop = FALSE], ols_fitfn, ols_predfn, alpha = alpha)
    cov[r] <- d$y[te] >= cf$lower && d$y[te] <= cf$upper
  }
  mean(cov)
}

test_that("couverture marginale >= 1-alpha (données gaussiennes)", {
  gen <- function(m) { X <- matrix(rnorm(m * 2), m, 2)
    list(X = X, y = 1 + X %*% c(2, -1) + rnorm(m)) }
  cov <- coverage_over_reps(gen, alpha = 0.1)
  expect_gte(cov, 0.90 - 0.02)                       # >= 1-alpha (à l'erreur MC près)
  expect_lte(cov, 0.90 + 0.05)
})

test_that("distribution-libre : couverture tient sous erreurs à queues lourdes", {
  gen <- function(m) { X <- matrix(rnorm(m * 2), m, 2)
    list(X = X, y = 1 + X %*% c(2, -1) + rt(m, df = 2)) }   # queues lourdes (t2)
  cov <- coverage_over_reps(gen, alpha = 0.1)
  expect_gte(cov, 0.90 - 0.02)                       # conforme reste valide
})

test_that("distribution-libre : couverture tient sous hétéroscédasticité", {
  gen <- function(m) { X <- matrix(rnorm(m * 2), m, 2)
    list(X = X, y = 1 + X %*% c(2, -1) + rnorm(m, sd = exp(0.6 * X[, 1]))) }
  cov <- coverage_over_reps(gen, alpha = 0.2)
  expect_gte(cov, 0.80 - 0.02)
})

test_that("un mauvais modèle donne des intervalles larges mais VALIDES", {
  set.seed(3)
  gen <- function(m) { X <- matrix(rnorm(m * 2), m, 2)
    list(X = X, y = sin(3 * X[, 1]) + X[, 2]^2 + rnorm(m, sd = 0.3)) }  # non linéaire
  # modèle linéaire (mal spécifié) : couverture toujours garantie
  cov <- coverage_over_reps(gen, alpha = 0.1, n = 300)
  expect_gte(cov, 0.90 - 0.02)
})
