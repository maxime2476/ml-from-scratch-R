# Tests de conformité — Module 1 (OLS + inférence). Référence : lm/summary.lm.
# Tolérance 1e-8.

make_df <- function(n = 100, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = runif(n),
                  g = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  d$y <- 1 - 1.5 * d$x1 + 0.5 * d$x2 + rnorm(n)
  d
}

test_that("ols_fit reproduit lm() : coefficients, fitted, résidus, leviers", {
  d <- make_df()
  fit <- ols_fit(y ~ x1 + x2 + x3, d)
  ref <- lm(y ~ x1 + x2 + x3, d)
  expect_equal(as.numeric(fit$coefficients), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(fit$fitted, as.numeric(fitted(ref)), tolerance = 1e-8)
  expect_equal(fit$residuals, as.numeric(residuals(ref)), tolerance = 1e-8)
  expect_equal(unname(fit$hat), as.numeric(hatvalues(ref)), tolerance = 1e-8)
  expect_equal(fit$df.residual, ref$df.residual)
})

test_that("vcov et sigma reproduisent lm()", {
  d <- make_df()
  fit <- ols_fit(y ~ x1 + x2 + x3, d)
  ref <- lm(y ~ x1 + x2 + x3, d)
  expect_equal(unname(fit$vcov), unname(vcov(ref)), tolerance = 1e-8)
  expect_equal(fit$sigma, summary(ref)$sigma, tolerance = 1e-8)
})

test_that("ols_summary reproduit summary.lm : se, t, p, R2, adj R2, F", {
  d <- make_df()
  sm <- ols_summary(ols_fit(y ~ x1 + x2 + x3, d))
  sref <- summary(lm(y ~ x1 + x2 + x3, d))
  expect_equal(sm$coefficients$se, unname(sref$coefficients[, "Std. Error"]), tolerance = 1e-8)
  expect_equal(sm$coefficients$t,  unname(sref$coefficients[, "t value"]),   tolerance = 1e-8)
  expect_equal(sm$coefficients$p_value, unname(sref$coefficients[, "Pr(>|t|)"]), tolerance = 1e-8)
  expect_equal(sm$r2, sref$r.squared, tolerance = 1e-8)
  expect_equal(sm$adj_r2, sref$adj.r.squared, tolerance = 1e-8)
  expect_equal(sm$fstatistic$value, unname(sref$fstatistic["value"]), tolerance = 1e-8)
  expect_equal(sm$fstatistic$p_value,
               pf(sref$fstatistic["value"], sref$fstatistic["numdf"],
                  sref$fstatistic["dendf"], lower.tail = FALSE)[[1]], tolerance = 1e-8)
})

test_that("gestion des facteurs (model.matrix) conforme à lm", {
  d <- make_df()
  fit <- ols_fit(y ~ x1 + g, d)
  ref <- lm(y ~ x1 + g, d)
  expect_equal(as.numeric(fit$coefficients), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(names(fit$coefficients), names(coef(ref)))
})

test_that("ols_confint reproduit confint.lm", {
  d <- make_df()
  fit <- ols_fit(y ~ x1 + x2 + x3, d)
  ref <- lm(y ~ x1 + x2 + x3, d)
  expect_equal(unname(ols_confint(fit, 0.95)), unname(confint(ref, level = 0.95)), tolerance = 1e-8)
  expect_equal(unname(ols_confint(fit, 0.90)), unname(confint(ref, level = 0.90)), tolerance = 1e-8)
})

test_that("ols_ftest : F d'une restriction = t^2 ; F conjoint = anova", {
  d <- make_df()
  fit <- ols_fit(y ~ x1 + x2 + x3, d)
  sm  <- ols_summary(fit)
  # restriction unique x3 = 0 -> F = t_x3^2
  R1 <- matrix(0, 1, 4); R1[1, 4] <- 1
  expect_equal(ols_ftest(fit, R1)$F, sm$coefficients$t[4]^2, tolerance = 1e-8)
  # restriction conjointe x2 = x3 = 0 -> comparer au F de comparaison de modèles
  R2 <- rbind(c(0, 0, 1, 0), c(0, 0, 0, 1))
  ft <- ols_ftest(fit, R2)
  ref_full <- lm(y ~ x1 + x2 + x3, d); ref_red <- lm(y ~ x1, d)
  Fref <- anova(ref_red, ref_full)$F[2]
  expect_equal(ft$F, Fref, tolerance = 1e-8)
})

test_that("Frisch-Waugh-Lovell : beta2 identique, vérification numérique", {
  d <- make_df()
  fit <- ols_fit(y ~ x1 + x2 + x3, d)
  X <- model.matrix(y ~ x1 + x2 + x3, d)
  # partialiser (Intercept, x1) ; bloc d'intérêt (x2, x3)
  b2 <- fwl_beta2(d$y, X[, c(1, 2)], X[, c(3, 4)])
  expect_equal(b2, as.numeric(fit$coefficients[c(3, 4)]), tolerance = 1e-8)
  # les résidus de la régression résidualisée = résidus de la régression complète
  resid_on_X1 <- function(z) z - X[, 1:2] %*% solve_ls_qr(X[, 1:2], z)$coefficients
  M1y <- resid_on_X1(d$y); M1X2 <- apply(X[, 3:4], 2, resid_on_X1)
  r_fwl <- as.numeric(M1y - M1X2 %*% b2)
  expect_equal(r_fwl, fit$residuals, tolerance = 1e-8)
})

test_that("R2 augmente mécaniquement quand on ajoute une variable de bruit", {
  set.seed(99)
  d <- make_df(120, seed = 99); d$noise <- rnorm(120)
  r2_small <- ols_summary(ols_fit(y ~ x1 + x2, d))$r2
  r2_big   <- ols_summary(ols_fit(y ~ x1 + x2 + noise, d))$r2
  expect_gte(r2_big, r2_small - 1e-12)   # non décroissant (Prop. 1.4)
})

test_that("sans constante : R2 non défini, coefficients corrects", {
  d <- make_df()
  fit <- ols_fit(y ~ x1 + x2 - 1, d)
  ref <- lm(y ~ x1 + x2 - 1, d)
  expect_equal(as.numeric(fit$coefficients), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_false(fit$has_intercept)
})
