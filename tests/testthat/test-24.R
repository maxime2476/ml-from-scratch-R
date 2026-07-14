# Tests — Module 24 (fonctions d'influence). Références : sandwich HC0 (M2),
# bootstrap (M17), GLM (M3). Unification vérifiée numériquement.

test_that("IC de l'OLS = sandwich HC0 (Prop. 24.1) à 1e-8", {
  set.seed(1); n <- 200
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n)); d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(n)
  X <- model.matrix(y ~ x1 + x2, d)
  inf <- influence_ols(X, d$y)
  hc0 <- vcov_hc(ols_fit(y ~ x1 + x2, d), "HC0")
  expect_lt(max(abs(inf$vcov - hc0)), 1e-8)
  expect_lt(max(abs(colMeans(inf$ic))), 1e-8)          # E[IC] = 0
})

test_that("variance IC = jackknife = bootstrap (pour la moyenne)", {
  set.seed(2); n <- 300; y <- rgamma(n, 2, 1)
  var_if <- var(y) / n
  jk <- jackknife(data.frame(y = y), function(dd) mean(dd$y))
  bt <- bootstrap(data.frame(y = y), function(dd) mean(dd$y), R = 3000, seed = 1)
  expect_equal(jk$var, var_if, tolerance = 1e-8)       # jackknife = IC (exact pour la moyenne)
  expect_lt(abs(var(bt$replicates) - var_if) / var_if, 0.1)  # bootstrap ~ IC
})

test_that("IC du MLE ~ variance du GLM (M3)", {
  set.seed(3); n <- 400; x <- rnorm(n); y <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm_irls(y ~ x, data.frame(y = y, x = x), family = "binomial")
  im <- influence_mle(fit)
  # sandwich robuste (IC) proche de la variance modèle (spéc. correcte)
  expect_lt(max(abs(sqrt(diag(im$vcov)) - fit$se)) / max(fit$se), 0.15)
})

test_that("estimateur en un pas : contraction vers le MLE (Newton)", {
  set.seed(4); n <- 500; x <- rnorm(n); y <- rpois(n, exp(0.7 * x))
  X <- cbind(1, x)
  mle <- coef(glm(y ~ x, family = poisson()))
  ic_at <- function(b) { mu <- exp(X %*% b); I <- crossprod(X * as.numeric(mu), X)
    (X * as.numeric(y - mu)) %*% t(solve(I)) * n }
  init <- mle + c(0.25, -0.25)                         # init dans un voisinage sqrt(n)
  os <- onestep(init, ic_at)
  # un pas de Newton contracte quadratiquement : bien plus proche du MLE
  expect_lt(max(abs(os - mle)), 0.3 * max(abs(init - mle)))
})

test_that("EIF de l'ATE : estimateur doublement robuste, couverture 95 %", {
  cover <- 0; R <- 300; z <- qnorm(0.975)
  for (r in seq_len(R)) {
    set.seed(1000 + r); n <- 500
    Z <- rnorm(n); e <- plogis(0.6 * Z); D <- rbinom(n, 1, e)
    mu0 <- 1 + Z; mu1 <- mu0 + 2                        # ATE = 2
    Y <- D * mu1 + (1 - D) * mu0 + rnorm(n)
    ea <- eif_ate(Y, D, mu1, mu0, e)
    cover <- cover + (abs(ea$ate - 2) <= z * ea$se)
  }
  expect_gt(cover / R, 0.92)
})
