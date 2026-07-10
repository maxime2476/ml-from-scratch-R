# Tests de conformité — Module 2 (hétéroscédasticité, robustesse, GLS/WLS).
# Références : sandwich::vcovHC / NeweyWest, lmtest::coeftest, lm(weights=),
# nlme::gls. Tolérance 1e-8.

skip_if_not_installed("sandwich")
skip_if_not_installed("lmtest")
skip_if_not_installed("nlme")

make_hetero <- function(n = 70, seed = 3) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$sig <- exp(0.5 * d$x1)                 # variance croissante avec x1
  d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(n, sd = d$sig)
  d
}

test_that("vcov_hc reproduit sandwich::vcovHC pour HC0-HC3", {
  d <- make_hetero()
  fit <- ols_fit(y ~ x1 + x2, d)
  ref <- lm(y ~ x1 + x2, d)
  for (ty in c("HC0", "HC1", "HC2", "HC3")) {
    expect_equal(unname(vcov_hc(fit, ty)),
                 unname(sandwich::vcovHC(ref, type = ty)), tolerance = 1e-8)
  }
})

test_that("coeftest_hc reproduit lmtest::coeftest (se, t, p) sous HC3", {
  d <- make_hetero()
  fit <- ols_fit(y ~ x1 + x2, d)
  ref <- lm(y ~ x1 + x2, d)
  ct  <- coeftest_hc(fit, vcov_hc(fit, "HC3"))
  ctr <- lmtest::coeftest(ref, vcov. = sandwich::vcovHC(ref, type = "HC3"))
  expect_equal(ct$se, unname(ctr[, "Std. Error"]), tolerance = 1e-8)
  expect_equal(ct$t,  unname(ctr[, "t value"]),   tolerance = 1e-8)
  expect_equal(ct$p_value, unname(ctr[, "Pr(>|t|)"]), tolerance = 1e-8)
})

test_that("wls_fit reproduit lm(weights=) : coefficients, vcov, sigma", {
  d <- make_hetero()
  w <- 1 / d$sig^2
  wf <- wls_fit(y ~ x1 + x2, d, weights = w)
  lw <- lm(y ~ x1 + x2, d, weights = w)
  expect_equal(as.numeric(wf$coefficients), as.numeric(coef(lw)), tolerance = 1e-8)
  expect_equal(unname(wf$vcov), unname(vcov(lw)), tolerance = 1e-8)
  expect_equal(wf$sigma, summary(lw)$sigma, tolerance = 1e-8)
  expect_equal(wf$residuals, as.numeric(residuals(lw)), tolerance = 1e-8)
})

test_that("gls_fit (Omega diagonale) = WLS et = nlme::gls varFixed", {
  d <- make_hetero()
  Omega <- diag(d$sig^2)
  gf <- gls_fit(y ~ x1 + x2, d, Omega = Omega)
  # équivalence GLS diagonal <-> WLS de poids 1/sig^2
  wf <- wls_fit(y ~ x1 + x2, d, weights = 1 / d$sig^2)
  expect_equal(as.numeric(gf$coefficients), as.numeric(wf$coefficients), tolerance = 1e-8)
  # vs nlme::gls avec structure de variance connue (varFixed)
  d$sig2 <- d$sig^2
  g_nlme <- nlme::gls(y ~ x1 + x2, data = d, weights = nlme::varFixed(~ sig2))
  expect_equal(as.numeric(gf$coefficients), as.numeric(coef(g_nlme)), tolerance = 1e-8)
})

test_that("gls_fit résout bien la forme fermée (X'Omega^-1 X)^-1 X'Omega^-1 y", {
  d <- make_hetero(50)
  set.seed(9)
  # Omega SPD non diagonale (corrélation AR(1))
  rho <- 0.4; n <- nrow(d)
  Omega <- rho^abs(outer(seq_len(n), seq_len(n), "-"))
  gf <- gls_fit(y ~ x1 + x2, d, Omega = Omega)
  X <- model.matrix(y ~ x1 + x2, d); y <- d$y
  Oinv <- solve(Omega)
  beta_ref <- as.numeric(solve(t(X) %*% Oinv %*% X, t(X) %*% Oinv %*% y))
  expect_equal(as.numeric(gf$coefficients), beta_ref, tolerance = 1e-8)
})

test_that("vcov_nw reproduit sandwich::NeweyWest (prewhite=FALSE, adjust=FALSE)", {
  d <- make_hetero()
  fit <- ols_fit(y ~ x1 + x2, d)
  ref <- lm(y ~ x1 + x2, d)
  for (L in c(1, 3, 5)) {
    expect_equal(unname(vcov_nw(fit, lag = L)),
                 unname(sandwich::NeweyWest(ref, lag = L, prewhite = FALSE,
                                            adjust = FALSE)), tolerance = 1e-8)
  }
})

test_that("Prop. 2.2 : E[e_i^2] = sigma^2 (1 - h_ii) sous homoscédasticité", {
  # Vérification numérique de la correction de levier par Monte Carlo court.
  set.seed(1)
  n <- 40
  X <- cbind(1, rnorm(n), rnorm(n))
  H <- X %*% solve(crossprod(X), t(X)); h <- diag(H)
  M <- diag(n) - H
  sigma <- 2
  R <- 20000
  e2 <- matrix(0, R, n)
  for (r in seq_len(R)) {
    eps <- rnorm(n, sd = sigma)
    e2[r, ] <- as.numeric(M %*% eps)^2
  }
  expect_equal(colMeans(e2), sigma^2 * (1 - h), tolerance = 0.05)
})
