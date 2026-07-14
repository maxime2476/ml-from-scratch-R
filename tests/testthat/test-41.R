# Tests — Module 41 (clustering avance). References : stats, dbscan.

ari <- function(a, b) {
  tab <- table(a, b); s <- sum(choose(tab, 2))
  a1 <- sum(choose(rowSums(tab), 2)); b1 <- sum(choose(colSums(tab), 2))
  n <- length(a); e <- a1 * b1 / choose(n, 2); (s - e) / ((a1 + b1) / 2 - e)
}
circles <- function(n = 200, seed = 1) {
  set.seed(seed); th <- runif(n, 0, 2 * pi); r <- c(rep(1, n / 2), rep(4, n / 2)) + rnorm(n, 0, 0.1)
  list(X = cbind(r * cos(th), r * sin(th)), y = rep(1:2, each = n / 2))
}

test_that("hierarchique : partition a k = cutree(hclust) (ARI = 1)", {
  set.seed(1); X <- rbind(matrix(rnorm(40), 20, 2), matrix(rnorm(40), 20, 2) + 4); D <- dist(X)
  for (m in c("complete", "single", "average"))
    expect_equal(ari(agglomerative(D, 3, m), cutree(hclust(D, m), 3)), 1, tolerance = 1e-8)
})

test_that("DBSCAN = dbscan::dbscan (meme partition)", {
  skip_if_not_installed("dbscan")
  set.seed(2); Y <- rbind(matrix(rnorm(60), 30, 2), matrix(rnorm(60), 30, 2) + 5)
  lh <- dbscan_fit(Y, eps = 1.0, minPts = 4); lr <- dbscan::dbscan(Y, eps = 1.0, minPts = 4)$cluster
  expect_true(all(outer(lh, lh, "==") == outer(lr, lr, "==")))
})

test_that("DBSCAN detecte le bruit (points isoles -> label 0)", {
  set.seed(3); Y <- rbind(matrix(rnorm(60), 30, 2), matrix(rnorm(60), 30, 2) + 6, c(20, 20), c(-15, 15))
  lab <- dbscan_fit(Y, eps = 1.0, minPts = 4)
  expect_equal(lab[61], 0L); expect_equal(lab[62], 0L)     # les 2 aberrants sont du bruit
})

test_that("spectral separe des cercles concentriques (la ou k-means echoue)", {
  d <- circles()
  sp <- spectral_clustering(d$X, 2, gamma = 1)
  km <- kmeans(d$X, 2, nstart = 10)$cluster
  expect_gt(ari(sp, d$y), 0.95)                            # spectral reussit
  expect_lt(ari(km, d$y), 0.2)                             # k-means echoue
})

test_that("agglomerative : nombre de groupes correct et couverture complete", {
  set.seed(4); X <- matrix(rnorm(60), 30, 2)
  a <- agglomerative(X, 4, "average")
  expect_equal(length(unique(a)), 4); expect_equal(length(a), 30)
})
