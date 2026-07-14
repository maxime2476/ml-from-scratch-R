# Tests — Module 47 (panel avance). References : construction + plm.

test_that("controle synthetique : recupere une combinaison convexe connue", {
  set.seed(1); Tt <- 30; J <- 8; pre <- 1:20
  Y0 <- matrix(rnorm(Tt * J), Tt, J); w_true <- c(0.5, 0.3, 0.2, rep(0, J - 3))
  Y1 <- as.numeric(Y0 %*% w_true)
  sc <- synthetic_control(Y1, Y0, pre)
  expect_lt(max(abs(sc$weights - w_true)), 1e-3)           # retrouve les poids
  expect_lt(sqrt(mean((Y1[pre] - sc$synthetic[pre])^2)), 1e-4)   # ajustement pre parfait
})

test_that("controle synthetique : poids sur le simplexe (>=0, somme 1)", {
  set.seed(2); Y0 <- matrix(rnorm(30 * 6), 30, 6); Y1 <- rnorm(30)
  sc <- synthetic_control(Y1, Y0, 1:20)
  expect_equal(sum(sc$weights), 1, tolerance = 1e-6)
  expect_true(all(sc$weights >= -1e-8))
})

test_that("controle synthetique detecte un effet post-traitement injecte", {
  set.seed(3); Tt <- 30; J <- 6; pre <- 1:20; post <- 21:30
  Y0 <- matrix(rnorm(Tt * J), Tt, J); w <- c(0.4, 0.35, 0.25, rep(0, J - 3))
  Y1 <- as.numeric(Y0 %*% w); Y1[post] <- Y1[post] + 2                 # effet de +2 en post
  sc <- synthetic_control(Y1, Y0, pre)
  expect_gt(mean(sc$effect), 1.5)                          # effet detecte (~ +2)
})

test_that("panel dynamique : l'IV corrige le biais de Nickell des effets fixes", {
  # sur plusieurs jeux, l'IV est ~ sans biais, le within nettement biaise vers le bas
  bias_fe <- bias_iv <- 0; R <- 40; rho <- 0.5; N <- 100; T <- 8
  for (r in 1:R) {
    set.seed(r); a <- rnorm(N); Y <- matrix(0, T, N)
    for (i in 1:N) { yi <- rnorm(1); for (t in 1:T) { yi <- rho * yi + a[i] + rnorm(1); Y[t, i] <- yi } }
    d <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N), y = as.numeric(Y))
    bias_fe <- bias_fe + (dynamic_panel_fe(d) - rho)
    bias_iv <- bias_iv + (dynamic_panel_iv(d)$rho - rho)
  }
  expect_lt(bias_fe / R, -0.1)                             # within : biais NEGATIF marque (Nickell)
  expect_lt(abs(bias_iv / R), 0.08)                        # IV : ~ sans biais
})

test_that("panel dynamique IV ~ plm::pgmm (Arellano-Bond)", {
  skip_if_not_installed("plm")
  set.seed(5); N <- 150; T <- 8; rho <- 0.6; a <- rnorm(N); Y <- matrix(0, T, N)
  for (i in 1:N) { yi <- rnorm(1); for (t in 1:T) { yi <- rho * yi + a[i] + rnorm(1); Y[t, i] <- yi } }
  d <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N), y = as.numeric(Y))
  iv <- dynamic_panel_iv(d)$rho
  pg <- tryCatch(as.numeric(coef(plm::pgmm(y ~ lag(y, 1) | lag(y, 2:99),
        data = plm::pdata.frame(d, index = c("id", "time")), effect = "individual", model = "onestep"))[1]),
        error = function(e) NA)
  skip_if(is.na(pg))
  expect_lt(abs(iv - pg), 0.15)                            # meme ordre de grandeur
})
