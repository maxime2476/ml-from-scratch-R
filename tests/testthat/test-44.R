# Tests — Module 44 (VAR). Reference : vars.

make_var <- function(n = 300, seed = 1) {
  set.seed(seed); Y <- matrix(0, n, 2)
  A1 <- matrix(c(0.5, 0.1, -0.2, 0.4), 2); A2 <- matrix(c(0.1, 0, 0.1, -0.1), 2)
  for (t in 3:n) Y[t, ] <- A1 %*% Y[t - 1, ] + A2 %*% Y[t - 2, ] + rnorm(2) * 0.5
  colnames(Y) <- c("y1", "y2"); Y
}

test_that("coefficients VAR = vars::VAR", {
  skip_if_not_installed("vars")
  Y <- make_var(); m <- var_fit(Y, 2); vr <- vars::VAR(Y, p = 2, type = "const")
  ch <- m$B[, 1]; ch_re <- c(ch[2], ch[3], ch[4], ch[5], ch[1])   # reordonner comme vars
  expect_lt(max(abs(ch_re - coef(vr)$y1[, 1])), 1e-8)
})

test_that("causalite de Granger = vars::causality", {
  skip_if_not_installed("vars")
  Y <- make_var(); m <- var_fit(Y, 2); vr <- vars::VAR(Y, p = 2, type = "const")
  g <- granger_test(m, cause = 2, effect = 1)
  expect_lt(abs(g$statistic - as.numeric(vars::causality(vr, cause = "y2")$Granger$statistic)), 1e-6)
})

test_that("IRF orthogonalisees = vars::irf", {
  skip_if_not_installed("vars")
  Y <- make_var(); m <- var_fit(Y, 2); vr <- vars::VAR(Y, p = 2, type = "const")
  ih <- var_irf(m, 10); ir <- vars::irf(vr, n.ahead = 10, ortho = TRUE, boot = FALSE)
  expect_lt(max(abs(ih[, , 1] - ir$irf$y1)), 1e-8)          # choc en y1
  expect_lt(max(abs(ih[, , 2] - ir$irf$y2)), 1e-8)          # choc en y2
})

test_that("Granger detecte la vraie direction (y2 -> y1, pas l'inverse fort)", {
  set.seed(2); n <- 500; Y <- matrix(0, n, 2)
  for (t in 2:n) { Y[t, 1] <- 0.3 * Y[t - 1, 1] + 0.6 * Y[t - 1, 2] + rnorm(1)   # y2 cause y1
                   Y[t, 2] <- 0.4 * Y[t - 1, 2] + rnorm(1) }                     # y1 ne cause pas y2
  m <- var_fit(Y, 1)
  expect_lt(granger_test(m, cause = 2, effect = 1)$p_value, 0.01)   # y2 -> y1 detecte
  expect_gt(granger_test(m, cause = 1, effect = 2)$p_value, 0.05)   # y1 -> y2 non detecte
})

test_that("IRF : le choc decroit (systeme stationnaire)", {
  Y <- make_var(); ih <- var_irf(var_fit(Y, 2), 20)
  # la reponse propre (y1 a un choc y1) decroit vers 0
  resp <- ih[, 1, 1]
  expect_lt(abs(resp[21]), abs(resp[1]))                    # decroissance
  expect_lt(abs(resp[21]), 0.1)
})
