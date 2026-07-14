# Tests — Module 42 (reduction de dimension). References : kernlab, fastICA.

ari <- function(a, b) { tab <- table(a, b); s <- sum(choose(tab, 2))
  a1 <- sum(choose(rowSums(tab), 2)); b1 <- sum(choose(colSums(tab), 2))
  n <- length(a); e <- a1 * b1 / choose(n, 2); (s - e) / ((a1 + b1) / 2 - e) }

test_that("kernel PCA : projection = kernlab::kpca (a un signe/echelle pres)", {
  skip_if_not_installed("kernlab")
  set.seed(1); X <- matrix(rnorm(60), 30, 2)
  kh <- kernel_pca(X, k = 2, gamma = 0.5)
  kr <- kernlab::kpca(X, kernel = "rbfdot", kpar = list(sigma = 0.5), features = 2)
  expect_gt(abs(cor(kh$proj[, 1], kernlab::rotated(kr)[, 1])), 0.999)
  expect_gt(abs(cor(kh$proj[, 2], kernlab::rotated(kr)[, 2])), 0.999)
})

test_that("FastICA recupere les sources independantes (cocktail party)", {
  set.seed(2); n <- 1000
  S <- cbind(sin((1:n) / 20), ((1:n) %% 50) / 50 - 0.5)     # 2 sources non gaussiennes
  X <- S %*% t(matrix(c(1, 0.6, 0.5, 1), 2))                 # melange
  ic <- ica_fastica(X, 2)
  # chaque source estimee correle a ~1 (en valeur absolue) avec une vraie source
  m <- abs(cor(ic$S, S))
  expect_gt(max(m[1, ]), 0.95); expect_gt(max(m[2, ]), 0.95)
})

test_that("NMF : reconstruction exacte d'un produit non negatif, W,H >= 0", {
  set.seed(3); Wt <- matrix(runif(20 * 3), 20, 3); Ht <- matrix(runif(3 * 15), 3, 15); V <- Wt %*% Ht
  nm <- nmf(V, 3, iter = 800)
  expect_lt(mean((V - nm$reconstruction)^2), 1e-4)          # reconstruit le produit
  expect_true(min(nm$W) >= 0 && min(nm$H) >= 0)             # non negativite
})

test_that("t-SNE : le plongement 2D preserve les clusters", {
  set.seed(4)
  X <- rbind(matrix(rnorm(40 * 8), 40, 8), matrix(rnorm(40 * 8), 40, 8) + 4,
             matrix(rnorm(40 * 8), 40, 8) - 4)
  y <- rep(1:3, each = 40)
  Y <- tsne(X, 2, perplexity = 30, iter = 300)
  expect_gt(ari(kmeans(Y, 3, nstart = 10)$cluster, y), 0.9)  # clusters recuperables en 2D
})

test_that("kernel PCA capte une structure non lineaire (cercles concentriques)", {
  set.seed(5); n <- 200; th <- runif(n, 0, 2 * pi); r <- c(rep(1, n/2), rep(4, n/2))
  X <- cbind(r * cos(th), r * sin(th)); y <- rep(1:2, each = n/2)
  kp <- kernel_pca(X, k = 2, gamma = 0.5)
  # la 1re composante kernel-PCA separe les rayons (les 2 cercles)
  expect_gt(abs(cor(kp$proj[, 1], y)), 0.9)
})
