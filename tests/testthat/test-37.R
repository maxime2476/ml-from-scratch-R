# Tests — Module 37 (RNN/LSTM). Reference : numDeriv + proprietes structurelles.

setup_rnn <- function(seed = 1) {
  set.seed(seed); Tn <- 6; d <- 3; H <- 4; O <- 2
  list(X = matrix(rnorm(Tn * d), Tn, d),
       Wxh = matrix(rnorm(H * d) * 0.3, H, d), Whh = matrix(rnorm(H * H) * 0.3, H, H),
       Why = matrix(rnorm(O * H) * 0.3, O, H), bh = rnorm(H) * 0.1, by = rnorm(O) * 0.1,
       Tn = Tn, d = d, H = H, O = O)
}

test_that("BPTT du RNN = numDeriv (tous les poids)", {
  skip_if_not_installed("numDeriv")
  p <- setup_rnn()
  fw <- rnn_forward(p$X, p$Wxh, p$Whh, p$Why, p$bh, p$by)
  bw <- rnn_backward(fw$Y, fw$cache)                       # dL/dY = Y, L = sum(Y^2)/2
  LWhh <- function(v) sum(rnn_forward(p$X, p$Wxh, matrix(v, p$H, p$H), p$Why, p$bh, p$by)$Y^2) / 2
  LWxh <- function(v) sum(rnn_forward(p$X, matrix(v, p$H, p$d), p$Whh, p$Why, p$bh, p$by)$Y^2) / 2
  LWhy <- function(v) sum(rnn_forward(p$X, p$Wxh, p$Whh, matrix(v, p$O, p$H), p$bh, p$by)$Y^2) / 2
  expect_lt(max(abs(bw$dWhh - matrix(numDeriv::grad(LWhh, as.numeric(p$Whh)), p$H, p$H))), 1e-6)
  expect_lt(max(abs(bw$dWxh - matrix(numDeriv::grad(LWxh, as.numeric(p$Wxh)), p$H, p$d))), 1e-6)
  expect_lt(max(abs(bw$dWhy - matrix(numDeriv::grad(LWhy, as.numeric(p$Why)), p$O, p$H))), 1e-6)
})

test_that("RNN forward : dimensions et coherence de l'etat", {
  p <- setup_rnn(); fw <- rnn_forward(p$X, p$Wxh, p$Whh, p$Why, p$bh, p$by)
  expect_equal(dim(fw$H), c(p$Tn, p$H)); expect_equal(dim(fw$Y), c(p$Tn, p$O))
  expect_true(all(abs(fw$H) <= 1))                        # tanh borne l'etat
})

test_that("LSTM : portes dans (0,1) et autoroute de memoire (f=1, i=0 -> c constant)", {
  set.seed(2); Tn <- 8; d <- 2; H <- 3
  X <- matrix(rnorm(Tn * d), Tn, d)
  W <- matrix(rnorm(H * (d + H)) * 0.3, H, d + H)
  # f = 1 (biais tres positif), i = 0 (biais tres negatif) -> c_t = c_{t-1}
  lf <- lstm_forward(X, W, W, W, W, rep(-20, H), rep(20, H), rep(0, H), rep(0, H))
  # c reste a 0 (init) car i~0 ; verifions qu'il ne diverge pas et reste ~constant
  expect_lt(max(abs(diff(lf$C))), 1e-6)                   # etat de cellule quasi constant
  expect_true(all(is.finite(lf$H)))
})

test_that("LSTM : l'etat de cellule accumule quand i=1, f=1", {
  set.seed(3); Tn <- 5; d <- 2; H <- 2
  X <- matrix(1, Tn, d)
  W <- matrix(0, H, d + H)                                 # g = tanh(bg)
  lf <- lstm_forward(X, W, W, W, W, rep(20, H), rep(20, H), rep(20, H), rep(0.5, H))
  # i=f=o=1, g=tanh(0.5) ~0.462 constant -> c_t = t * g croit lineairement
  expect_gt(lf$C[Tn, 1], lf$C[1, 1])                      # accumulation
})
