# Tests — Module 17 (bootstrap). Référence : boot::boot / boot.ci.

skip_if_not_installed("boot")

test_that("Prop. 17.1 : se bootstrap de la moyenne ~ sd(x)/sqrt(n)", {
  set.seed(1); x <- rnorm(80, 5, 2)
  bt <- bootstrap(x, mean, R = 5000, seed = 1)
  expect_equal(bt$se, sd(x) / sqrt(length(x)), tolerance = 0.06)   # (n-1)/n près
  expect_equal(bt$t0, mean(x))
})

test_that("percentile = quantiles des répliques (définitionnel)", {
  set.seed(2); x <- rgamma(60, 2, 1)
  bt <- bootstrap(x, mean, R = 3000, seed = 2)
  expect_equal(boot_ci(bt, type = "percentile"),
               unname(quantile(bt$replicates, c(.025, .975), type = 6)))
})

test_that("boot_ci reproduit boot::boot.ci (percentile, basique, BCa)", {
  set.seed(3); x <- rgamma(70, 2, 1)                 # asymétrique
  statfn <- function(d, i) mean(d[i])
  bo <- boot::boot(x, statfn, R = 6000)
  ref <- boot::boot.ci(bo, type = c("perc", "basic", "bca"))
  bt <- bootstrap(x, mean, R = 6000, seed = 10)
  w <- diff(ref$percent[4:5])                         # largeur ~ échelle de tolérance
  expect_lt(max(abs(boot_ci(bt, type = "percentile") - ref$percent[4:5])), 0.10 * w)
  expect_lt(max(abs(boot_ci(bt, type = "basic") - ref$basic[4:5])), 0.10 * w)
  expect_lt(max(abs(boot_ci(bt, type = "bca") - ref$bca[4:5])), 0.15 * w)
})

test_that("boot_lm : résidus ~ se classique ; pairs ~ se robuste (HC)", {
  set.seed(4); n <- 200
  d <- data.frame(x = rnorm(n)); d$y <- 1 + 2 * d$x + rnorm(n, sd = exp(0.5 * d$x))  # hétéro
  se_classic <- summary(lm(y ~ x, d))$coefficients[, "Std. Error"]
  se_hc <- sqrt(diag(sandwich::vcovHC(lm(y ~ x, d), type = "HC0")))
  br <- boot_lm(y ~ x, d, R = 3000, method = "residual", seed = 1)
  bp <- boot_lm(y ~ x, d, R = 3000, method = "pairs", seed = 1)
  # résidus : proche de la se classique (i.i.d.) pour la pente
  expect_lt(abs(br$se["x"] - se_classic["x"]) / se_classic["x"], 0.15)
  # pairs : proche de la se robuste HC (capte l'hétéroscédasticité)
  expect_lt(abs(bp$se["x"] - se_hc["x"]) / se_hc["x"], 0.15)
})
