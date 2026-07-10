# Tests — Module 7 (KNN, fléau de la dimension). Référence : class::knn.

skip_if_not_installed("class")

test_that("knn_classify reproduit class::knn (k impair, 2 classes, continu)", {
  set.seed(3)
  ntr <- 250; nte <- 120; p <- 4
  Xtr <- matrix(rnorm(ntr * p), ntr, p); Xte <- matrix(rnorm(nte * p), nte, p)
  ytr <- factor(ifelse(Xtr[, 1] + Xtr[, 2] + rnorm(ntr) > 0, "pos", "neg"))
  for (k in c(1, 3, 5, 7, 15)) {
    mine <- knn_classify(Xtr, ytr, Xte, k)
    ref  <- class::knn(Xtr, Xte, ytr, k = k)
    expect_identical(as.character(mine), as.character(ref))
  }
})

test_that("knn_regression : k=1 -> voisin le plus proche ; k=n -> moyenne globale", {
  set.seed(1)
  Xtr <- matrix(c(0, 1, 2, 10), 4, 1); ytr <- c(5, 6, 7, 100)
  # test point 1.1 : voisin le plus proche = ligne 2 (x=1) -> y=6
  expect_equal(knn_regression(Xtr, ytr, matrix(1.1, 1, 1), k = 1), 6)
  # k = n : moyenne de tout y quel que soit x
  expect_equal(knn_regression(Xtr, ytr, matrix(3, 1, 1), k = 4), mean(ytr))
})

test_that("knn_regression : moyenne locale correcte sur exemple contrôlé", {
  Xtr <- matrix(c(0, 1, 2, 3, 4), 5, 1); ytr <- c(0, 2, 4, 6, 8)
  # point 1.4 : 3 voisins les plus proches = x in {1,2,0} -> y in {2,4,0}, moyenne 2
  expect_equal(knn_regression(Xtr, ytr, matrix(1.4, 1, 1), k = 3), mean(c(2, 4, 0)))
})

test_that("edge_length : formule r^(1/p)", {
  expect_equal(edge_length(0.01, 10), 0.01^(1 / 10))
  expect_equal(edge_length(1, 5), 1)               # toute la donnée -> arête 1
  expect_gt(edge_length(0.01, 100), edge_length(0.01, 10))  # croît avec p
})

test_that("concentration des distances : cv(d^2) décroît en ~1/sqrt(p)", {
  cv2  <- distance_concentration(300, 2,  seed = 1)$cv_d2
  cv50 <- distance_concentration(300, 50, seed = 1)$cv_d2
  cv200 <- distance_concentration(300, 200, seed = 1)$cv_d2
  expect_gt(cv2, cv50)
  expect_gt(cv50, cv200)
  # scaling ~ 1/sqrt(p) : cv2/cv200 ~ sqrt(200/2) = 10 (tolérance large)
  expect_lt(abs(cv2 / cv200 / sqrt(200 / 2) - 1), 0.3)
})

test_that("concentration : le contraste (Dmax-Dmin)/Dmin s'effondre avec p", {
  c2   <- distance_concentration(300, 2,   seed = 2)$contrast
  c100 <- distance_concentration(300, 100, seed = 2)$contrast
  expect_gt(c2, c100)
  expect_lt(c100, 1)                                # voisins quasi équidistants
})
