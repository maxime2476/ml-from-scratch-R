# Tests — Module 32 (non parametrique + RDD). References : stats, rdrobust.

test_that("KDE = stats::density (noyau gaussien)", {
  set.seed(1); x <- rnorm(300); g <- seq(-3, 3, length.out = 40); h <- 0.4
  dref <- density(x, bw = h, kernel = "gaussian", n = 1024)
  di <- approx(dref$x, dref$y, g)$y
  expect_equal(kde(x, g, h), di, tolerance = 1e-3)          # density : binning FFT
})

test_that("Nadaraya-Watson = stats::ksmooth (avec conversion de fenetre)", {
  set.seed(2); x <- sort(runif(300, -3, 3)); y <- sin(x) + 0.3 * rnorm(300)
  bw <- 1.0; h <- 0.25 * bw / qnorm(0.75); x0 <- seq(-2, 2, length.out = 20)
  ks <- ksmooth(x, y, "normal", bandwidth = bw, x.points = x0)
  expect_equal(nadaraya_watson(x, y, x0, h), ks$y, tolerance = 1e-3)
})

test_that("regression locale = WLS ponderee (definition) et recupere la fonction", {
  set.seed(3); x <- sort(runif(400, -3, 3)); y <- sin(x) + 0.2 * rnorm(400)
  # coincide avec la definition (WLS locale)
  g <- 0.5; w <- dnorm((g - x) / 0.5); X <- cbind(1, x - g)
  b0 <- solve(crossprod(X * w, X), crossprod(X * w, y))[1]
  expect_equal(local_linear(x, y, g, 0.5), as.numeric(b0), tolerance = 1e-10)
  # recupere sin sur l'interieur
  x0 <- seq(-2, 2, length.out = 30)
  expect_lt(sqrt(mean((local_linear(x, y, x0, 0.4) - sin(x0))^2)), 0.1)
})

test_that("bw_loocv choisit une fenetre raisonnable", {
  set.seed(4); x <- sort(runif(200, -3, 3)); y <- sin(x) + 0.3 * rnorm(200)
  b <- bw_loocv(x, y, seq(0.1, 1.5, by = 0.1))
  expect_true(b$bw > 0.1 && b$bw < 1.5)
  expect_equal(length(b$cv), 15)
})

test_that("RDD : recupere le saut vrai et ~ rdrobust", {
  set.seed(5); x <- runif(600, -1, 1); tau <- 1.5
  y <- 0.5 * x + tau * (x >= 0) + 0.2 * rnorm(600)
  r <- rdd(y, x, cutoff = 0, bw = 0.3)
  expect_lt(abs(r$tau - tau), 0.15)                         # recupere le saut
  skip_if_not_installed("rdrobust")
  rr <- rdrobust::rdrobust(y, x, c = 0)
  rh <- rdd(y, x, cutoff = 0, bw = rr$bws[1, 1])            # meme fenetre, noyau tri
  expect_lt(abs(rh$tau - rr$Estimate[1, "tau.us"]), 0.1)    # ~ rdrobust
})
