# Tests — Module 11 (non supervisé). Références : prcomp, stats::kmeans, mclust.

skip_if_not_installed("mclust")
skip_if_not_installed("MASS")

test_that("pca_fit reproduit prcomp (sdev, rotation au signe près, var expliquée)", {
  set.seed(1)
  X <- matrix(rnorm(200 * 4), 200, 4); X[, 2] <- X[, 1] + 0.3 * rnorm(200)
  pc <- pca_fit(X); pr <- prcomp(X)
  expect_equal(pc$sdev, pr$sdev, tolerance = 1e-8)
  expect_equal(abs(pc$rotation), abs(unname(pr$rotation)), tolerance = 1e-8)  # signe libre
  expect_equal(pc$var_explained, pr$sdev^2 / sum(pr$sdev^2), tolerance = 1e-8)
})

test_that("ACP : les deux voies coïncident (Prop. 11.2)", {
  set.seed(2)
  X <- matrix(rnorm(150 * 3), 150, 3)
  pc <- pca_fit(X)
  eig <- eigen(cov(X))                                # voie 1 : vecteurs propres de S
  expect_equal(pc$sdev^2, eig$values, tolerance = 1e-8)              # lambda_j = d_j^2/(n-1)
  expect_equal(abs(pc$rotation), abs(eig$vectors), tolerance = 1e-8) # mêmes axes
})

test_that("kmeans_fit reproduit stats::kmeans (init identique, Lloyd)", {
  set.seed(2)
  X <- rbind(matrix(rnorm(100, 0), 50, 2), matrix(rnorm(100, 4), 50, 2))
  init <- X[c(5, 80), ]
  mine <- kmeans_fit(X, 2, centers = init)
  ref  <- kmeans(X, centers = init, algorithm = "Lloyd", iter.max = 100)
  expect_equal(mine$tot_withinss, ref$tot.withinss, tolerance = 1e-8)
  expect_identical(mine$cluster, ref$cluster)
})

test_that("k-means : l'inertie décroît de façon monotone (Prop. 11.3)", {
  set.seed(5)
  X <- rbind(matrix(rnorm(60, 0), 30, 2), matrix(rnorm(60, 5), 30, 2))
  # inertie après chaque itération, via kmeans à iter.max croissant
  inerties <- sapply(1:6, function(it) {
    kmeans_fit(X, 2, centers = X[c(1, 60), ], max_iter = it)$tot_withinss
  })
  expect_true(all(diff(inerties) <= 1e-8))
})

test_that("gmm_loglik évalué aux paramètres de mclust = loglik de mclust", {
  set.seed(3)
  Xm <- rbind(MASS::mvrnorm(80, c(0, 0), diag(2)),
              MASS::mvrnorm(80, c(5, 5), matrix(c(1, 0.5, 0.5, 1), 2)),
              MASS::mvrnorm(80, c(0, 6), diag(2)))
  mc <- mclust::Mclust(Xm, G = 3, modelNames = "VVV", verbose = FALSE)
  piM <- mc$parameters$pro
  muM <- lapply(1:3, function(k) mc$parameters$mean[, k])
  SigM <- lapply(1:3, function(k) mc$parameters$variance$sigma[, , k])
  expect_equal(gmm_loglik(Xm, piM, muM, SigM), mc$loglik, tolerance = 1e-6)
})

test_that("em_gmm atteint (au moins) la log-vraisemblance de mclust sur données séparées", {
  set.seed(3)
  Xm <- rbind(MASS::mvrnorm(80, c(0, 0), diag(2)),
              MASS::mvrnorm(80, c(6, 6), diag(2)),
              MASS::mvrnorm(80, c(0, 7), diag(2)))
  mc <- mclust::Mclust(Xm, G = 3, modelNames = "VVV", verbose = FALSE)
  em <- em_gmm(Xm, 3, seed = 1)
  # responsabilités : lignes sommant à 1 ; poids sommant à 1
  expect_equal(rowSums(em$gamma), rep(1, nrow(Xm)), tolerance = 1e-8)
  expect_equal(sum(em$pi), 1, tolerance = 1e-8)
  # EM converge vers un optimum au moins aussi bon (à tolérance numérique près)
  expect_gt(em$loglik, mc$loglik - 1)
})

test_that("k-means = limite EM : responsabilités quasi dures si clusters séparés", {
  set.seed(7)
  Xm <- rbind(MASS::mvrnorm(60, c(0, 0), 0.3 * diag(2)),
              MASS::mvrnorm(60, c(8, 8), 0.3 * diag(2)))
  em <- em_gmm(Xm, 2, seed = 1)
  expect_gt(mean(apply(em$gamma, 1, max)), 0.99)     # affectations quasi dures
})
