# Tests — Module 21 (panel / effets fixes). Références : lm (LSDV), plm, vcovCL.

skip_if_not_installed("plm")

make_panel <- function(N = 100, T = 6, seed = 5) {
  set.seed(seed)
  id <- rep(seq_len(N), each = T); tt <- rep(seq_len(T), N)
  alpha <- rnorm(N)
  x <- 0.7 * alpha[id] + rnorm(N * T)                 # x corrélé à l'effet fixe
  x2 <- rnorm(N * T)
  u <- rep(rnorm(N), each = T) * 0.5 + rnorm(N * T)   # corrélation intra-unité
  y <- 1 + 2 * x - 1 * x2 + alpha[id] + u
  data.frame(id = id, t = tt, x = x, x2 = x2, y = y)
}

test_that("within = LSDV (lm avec indicatrices) — Prop. 21.1", {
  d <- make_panel()
  fe <- fe_within(y ~ x + x2, d, "id")
  lsdv <- coef(lm(y ~ x + x2 + factor(id), d))[c("x", "x2")]
  expect_equal(as.numeric(fe$coefficients), as.numeric(lsdv), tolerance = 1e-8)
})

test_that("within = plm(model = 'within')", {
  d <- make_panel()
  fe <- fe_within(y ~ x + x2, d, "id")
  pm <- plm::plm(y ~ x + x2, data = d, index = c("id", "t"), model = "within")
  expect_equal(as.numeric(fe$coefficients), as.numeric(coef(pm)[c("x", "x2")]), tolerance = 1e-8)
})

test_that("SE groupées ~ plm (arellano, HC1)", {
  skip_if_not_installed("sandwich")
  d <- make_panel()
  fe <- fe_within(y ~ x + x2, d, "id")
  pm <- plm::plm(y ~ x + x2, data = d, index = c("id", "t"), model = "within")
  se_plm <- sqrt(diag(plm::vcovHC(pm, method = "arellano", type = "HC1")))
  expect_lt(max(abs(fe$se_cluster[c("x", "x2")] - se_plm[c("x", "x2")]) / se_plm[c("x", "x2")]), 0.05)
})

test_that("les SE groupées dépassent les SE naïves sous corrélation SÉRIELLE", {
  # x ET epsilon en AR(1) au sein de l'unité -> le clustering inflate la SE.
  set.seed(2); N <- 100; T <- 10; rho <- 0.8
  ar1 <- function() { e <- numeric(T); e[1] <- rnorm(1)
    for (t in 2:T) e[t] <- rho * e[t - 1] + rnorm(1); e }
  id <- rep(seq_len(N), each = T)
  x <- as.numeric(sapply(seq_len(N), function(i) ar1()))
  eps <- as.numeric(sapply(seq_len(N), function(i) ar1()))
  y <- 1 + 2 * x + rnorm(N)[id] + eps
  d <- data.frame(id = id, x = x, y = y)
  fe <- fe_within(y ~ x, d, "id")
  expect_gt(fe$se_cluster["x"], fe$se["x"])          # clustered > naïve
})

test_that("l'OLS groupé est biaisé (alpha corrélé), les effets fixes corrigent", {
  d <- make_panel(N = 300)
  b_pool <- coef(lm(y ~ x + x2, d))["x"]              # ignore alpha
  b_fe <- fe_within(y ~ x + x2, d, "id")$coefficients["x"]
  expect_lt(abs(b_fe - 2), abs(b_pool - 2))           # FE plus proche du vrai
})
