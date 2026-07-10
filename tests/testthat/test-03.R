# Tests de conformité — Module 3 (GLM par IRLS). Référence : glm(), anova.glm.
# Tolérance 1e-8.

make_glm_data <- function(n = 250, seed = 5) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = runif(n))
  eta <- -0.5 + 1.2 * d$x1 - 0.8 * d$x2 + 0.5 * d$x3
  d$yb <- rbinom(n, 1, plogis(eta))
  d$yp <- rpois(n, exp(0.3 + 0.5 * d$x1 - 0.4 * d$x2))
  d
}

test_that("logistique : glm_irls reproduit glm(binomial)", {
  d <- make_glm_data()
  fit <- glm_irls(yb ~ x1 + x2 + x3, d, "binomial")
  ref <- glm(yb ~ x1 + x2 + x3, d, family = binomial)
  expect_equal(as.numeric(fit$coefficients), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(unname(fit$se), unname(summary(ref)$coefficients[, "Std. Error"]), tolerance = 1e-8)
  expect_equal(unname(fit$vcov), unname(vcov(ref)), tolerance = 1e-8)
  expect_equal(fit$deviance, ref$deviance, tolerance = 1e-8)
  expect_equal(fit$null.deviance, ref$null.deviance, tolerance = 1e-8)
  expect_equal(unname(fit$fitted), unname(fitted(ref)), tolerance = 1e-8)
  expect_equal(fit$iter, ref$iter)
})

test_that("Poisson : glm_irls reproduit glm(poisson)", {
  d <- make_glm_data()
  fit <- glm_irls(yp ~ x1 + x2, d, "poisson")
  ref <- glm(yp ~ x1 + x2, d, family = poisson)
  expect_equal(as.numeric(fit$coefficients), as.numeric(coef(ref)), tolerance = 1e-8)
  expect_equal(unname(fit$se), unname(summary(ref)$coefficients[, "Std. Error"]), tolerance = 1e-8)
  expect_equal(fit$deviance, ref$deviance, tolerance = 1e-8)
  expect_equal(fit$null.deviance, ref$null.deviance, tolerance = 1e-8)
  expect_equal(unname(fit$fitted), unname(fitted(ref)), tolerance = 1e-8)
  expect_equal(fit$iter, ref$iter)
})

test_that("loglik cohérent avec logLik(glm)", {
  d <- make_glm_data()
  fb <- glm_irls(yb ~ x1 + x2, d, "binomial")
  rb <- glm(yb ~ x1 + x2, d, family = binomial)
  expect_equal(fb$loglik, as.numeric(logLik(rb)), tolerance = 1e-8)
  fp <- glm_irls(yp ~ x1 + x2, d, "poisson")
  rp <- glm(yp ~ x1 + x2, d, family = poisson)
  expect_equal(fp$loglik, as.numeric(logLik(rp)), tolerance = 1e-8)
})

test_that("test de Wald = z^2 de summary(glm) pour une restriction", {
  d <- make_glm_data()
  fit <- glm_irls(yb ~ x1 + x2 + x3, d, "binomial")
  ref <- glm(yb ~ x1 + x2 + x3, d, family = binomial)
  R <- matrix(c(0, 0, 1, 0), 1, 4)   # x2 = 0
  z2 <- (summary(ref)$coefficients["x2", "z value"])^2
  expect_equal(wald_test(fit, R)$statistic, z2, tolerance = 1e-8)
})

test_that("test LR reproduit anova(glm, test='LRT')", {
  d <- make_glm_data()
  full <- glm_irls(yb ~ x1 + x2 + x3, d, "binomial")
  red  <- glm_irls(yb ~ x1, d, "binomial")
  lr <- lr_test(full, red)
  a <- anova(glm(yb ~ x1, d, family = binomial),
             glm(yb ~ x1 + x2 + x3, d, family = binomial), test = "LRT")
  expect_equal(lr$statistic, a$Deviance[2], tolerance = 1e-8)
  expect_equal(lr$df, a$Df[2])
  expect_equal(lr$p_value, a$`Pr(>Chi)`[2], tolerance = 1e-8)
})

test_that("test du score reproduit anova(glm, test='Rao')", {
  d <- make_glm_data()
  full <- glm_irls(yb ~ x1 + x2 + x3, d, "binomial")
  red  <- glm_irls(yb ~ x1, d, "binomial")
  sc <- score_test(full, red)
  a <- anova(glm(yb ~ x1, d, family = binomial),
             glm(yb ~ x1 + x2 + x3, d, family = binomial), test = "Rao")
  expect_equal(sc$statistic, a$Rao[2], tolerance = 1e-8)
  expect_equal(sc$df, a$Df[2])
})

test_that("Poisson : LR et score reproduisent anova.glm", {
  d <- make_glm_data()
  full <- glm_irls(yp ~ x1 + x2, d, "poisson")
  red  <- glm_irls(yp ~ x1, d, "poisson")
  aL <- anova(glm(yp ~ x1, d, family = poisson),
              glm(yp ~ x1 + x2, d, family = poisson), test = "LRT")
  aR <- anova(glm(yp ~ x1, d, family = poisson),
              glm(yp ~ x1 + x2, d, family = poisson), test = "Rao")
  expect_equal(lr_test(full, red)$statistic, aL$Deviance[2], tolerance = 1e-8)
  expect_equal(score_test(full, red)$statistic, aR$Rao[2], tolerance = 1e-8)
})
