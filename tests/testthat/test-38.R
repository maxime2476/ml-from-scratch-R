# Tests — Module 38 (attention/Transformer). Reference : definition + proprietes.

test_that("softmax : lignes positives sommant a 1", {
  set.seed(1); X <- matrix(rnorm(20), 4, 5)
  S <- softmax_rows(X)
  expect_equal(rowSums(S), rep(1, 4), tolerance = 1e-12)
  expect_true(all(S > 0))
})

test_that("attention = definition ; poids sommant a 1", {
  set.seed(2); Tn <- 5; dk <- 4; dv <- 3
  Q <- matrix(rnorm(Tn * dk), Tn, dk); K <- matrix(rnorm(Tn * dk), Tn, dk); V <- matrix(rnorm(Tn * dv), Tn, dv)
  a <- attention(Q, K, V)
  W <- softmax_rows((Q %*% t(K)) / sqrt(dk))
  expect_equal(a$weights, W, tolerance = 1e-12)
  expect_equal(a$out, W %*% V, tolerance = 1e-12)
  expect_equal(rowSums(a$weights), rep(1, Tn), tolerance = 1e-12)
})

test_that("masque causal : la position i ne voit pas j > i", {
  set.seed(3); Tn <- 6; d <- 4
  Q <- matrix(rnorm(Tn * d), Tn, d); K <- matrix(rnorm(Tn * d), Tn, d); V <- matrix(rnorm(Tn * d), Tn, d)
  a <- attention(Q, K, V, mask = TRUE)
  expect_lt(max(abs(a$weights[upper.tri(a$weights)])), 1e-12)
  expect_equal(rowSums(a$weights), rep(1, Tn), tolerance = 1e-12)
})

test_that("self-attention permutation-equivariante ; l'encodage positionnel la brise", {
  set.seed(4); Tn <- 5; d <- 4
  X <- matrix(rnorm(Tn * d), Tn, d); perm <- c(3, 1, 5, 2, 4)
  sa <- function(Z) attention(Z, Z, Z)$out                 # self-attention Q=K=V=Z
  # sans encodage : permuter TOUTE la sequence permute identiquement les sorties
  expect_lt(max(abs(sa(X[perm, ]) - sa(X)[perm, ])), 1e-10)
  # avec encodage positionnel (lie a la POSITION), l'equivariance est brisee
  PE <- positional_encoding(Tn, d)
  expect_gt(max(abs(sa(X[perm, ] + PE) - sa(X + PE)[perm, ])), 1e-6)
})

test_that("attention comme recherche associative : une requete recupere sa valeur", {
  d <- 8; K <- diag(d)[1:4, ]; V <- matrix(1:4, 4, 1) * 10   # 4 cles orthonormees, valeurs distinctes
  q <- K[2, ] * 20                                           # requete alignee sur la cle 2 (forte)
  a <- attention(matrix(q, 1), K, V)
  expect_gt(a$weights[1, 2], 0.99)                          # poids concentre sur la cle 2
  expect_lt(abs(a$out[1, 1] - 20), 1)                       # recupere ~ V_2 = 20
})

test_that("multi-tetes et positional encoding : dimensions et formule", {
  set.seed(5); Tn <- 6; d <- 8; X <- matrix(rnorm(Tn * d), Tn, d)
  W <- matrix(rnorm(d * d) * 0.1, d, d)
  mh <- multi_head_attention(X, W, W, W, diag(d), n_heads = 4)
  expect_equal(dim(mh$out), c(Tn, d)); expect_equal(length(mh$weights), 4)
  PE <- positional_encoding(10, 6)
  expect_equal(PE[3, 1], sin(2), tolerance = 1e-12)         # pos=2, i=0 : sin(2)
  expect_equal(PE[3, 2], cos(2), tolerance = 1e-12)
})
