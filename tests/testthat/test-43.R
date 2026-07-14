# Tests — Module 43 (bandits/RL). Reference : theorie (regret sous-lineaire) + value iteration.

test_that("UCB et Thompson : regret sous-lineaire, bien inferieur au hasard", {
  means <- c(0.3, 0.5, 0.7, 0.55, 0.4); T <- 5000
  ru <- bandit_ucb(means, T, seed = 1); rt <- bandit_thompson(means, T, seed = 1)
  rand <- (max(means) - mean(means)) * T
  expect_lt(ru$regret[T], 0.4 * rand)        # bien mieux que le hasard
  expect_lt(rt$regret[T], 0.4 * rand)
  expect_lt(ru$regret[T] / T, 0.1)           # regret/T -> 0 (sous-lineaire)
  expect_lt(rt$regret[T] / T, 0.1)
})

test_that("les bandits convergent vers le meilleur bras", {
  means <- c(0.2, 0.4, 0.8, 0.5); T <- 4000
  ru <- bandit_ucb(means, T, seed = 2); rt <- bandit_thompson(means, T, seed = 2)
  # la seconde moitie du temps est majoritairement passee sur le bras 3
  expect_gt(mean(ru$arms[(T/2):T] == 3), 0.8)
  expect_gt(mean(rt$arms[(T/2):T] == 3), 0.9)
})

test_that("le regret croit sous-lineairement (pente log-log < 1)", {
  means <- c(0.3, 0.6, 0.5); ru <- bandit_ucb(means, 8000, seed = 3)
  ts <- c(500, 1000, 2000, 4000, 8000); reg <- ru$regret[ts]
  slope <- coef(lm(log(reg) ~ log(ts)))[2]
  expect_lt(slope, 1)                        # sous-lineaire (theorie : ~ log T)
})

test_that("value iteration : point fixe de Bellman", {
  set.seed(4); S <- 4; A <- 2; P <- array(0, c(S, A, S))
  for (s in 1:S) for (a in 1:A) { p <- runif(S); P[s, a, ] <- p / sum(p) }
  R <- matrix(runif(S * A), S, A); vi <- value_iteration(P, R, 0.9)
  # verifier l'equation de Bellman
  V <- vi$V; Qcheck <- R + 0.9 * apply(P, c(1, 2), function(p) sum(p * V))
  expect_lt(max(abs(vi$Q - Qcheck)), 1e-6)
})

test_that("Q-learning converge vers la politique optimale (= value iteration)", {
  set.seed(5); S <- 4; A <- 3; P <- array(0, c(S, A, S))
  for (s in 1:S) for (a in 1:A) { p <- runif(S); P[s, a, ] <- p / sum(p) }
  best_a <- c(2, 1, 3, 2)                     # action optimale CLAIRE par etat
  R <- matrix(0.1, S, A); for (s in 1:S) R[s, best_a[s]] <- 1   # ecart net de recompense
  vi <- value_iteration(P, R, 0.9)
  ql <- q_learning(P, R, 0.9, episodes = 20000, alpha = 0.1, epsilon = 0.2, seed = 6)
  expect_equal(ql$policy, vi$policy)         # meme politique optimale
  expect_lt(max(abs(ql$Q - vi$Q)), 0.35)     # Q proches (alpha constant : residu)
})
