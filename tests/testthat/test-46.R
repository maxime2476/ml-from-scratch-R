# Tests — Module 46 (survie). Reference : survival.

make_surv <- function(n = 200, beta = 0.8, seed = 1) {
  set.seed(seed); x <- rnorm(n); tt <- rexp(n, rate = exp(beta * x) * 0.1)
  cens <- runif(n, 0, 15); list(time = pmin(tt, cens), event = as.integer(tt <= cens), x = x)
}

test_that("Kaplan-Meier = survival::survfit", {
  skip_if_not_installed("survival")
  d <- make_surv(); km <- kaplan_meier(d$time, d$event)
  sf <- survival::survfit(survival::Surv(d$time, d$event) ~ 1)
  sr <- summary(sf, times = km$time)$surv
  expect_lt(max(abs(km$surv - sr)), 1e-8)
})

test_that("Cox : coefficient = survival::coxph (Breslow)", {
  skip_if_not_installed("survival")
  d <- make_surv()
  ch <- cox_ph(d$time, d$event, d$x)
  cr <- survival::coxph(survival::Surv(d$time, d$event) ~ d$x, ties = "breslow")
  expect_lt(abs(as.numeric(ch$coefficients) - as.numeric(coef(cr))), 1e-3)
  expect_lt(abs(as.numeric(ch$se) - sqrt(diag(vcov(cr)))), 1e-3)
})

test_that("log-rank = survival::survdiff", {
  skip_if_not_installed("survival")
  d <- make_surv(); grp <- as.integer(d$x > 0)
  lh <- logrank_test(d$time, d$event, grp)
  ls <- survival::survdiff(survival::Surv(d$time, d$event) ~ grp)
  expect_lt(abs(lh$statistic - ls$chisq), 1e-6)
})

test_that("Cox recupere le vrai rapport de risque", {
  d <- make_surv(n = 1000, beta = 1.2, seed = 3)
  ch <- cox_ph(d$time, d$event, d$x)
  expect_lt(abs(as.numeric(ch$coefficients) - 1.2), 0.15)   # ~ vrai beta
})

test_that("Kaplan-Meier : survie decroissante dans [0,1], gere la censure", {
  d <- make_surv(); km <- kaplan_meier(d$time, d$event)
  expect_true(all(km$surv >= 0 & km$surv <= 1))
  expect_true(all(diff(km$surv) <= 1e-12))                 # decroissante
})
