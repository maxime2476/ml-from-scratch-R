# =============================================================================
# Module 43 — Bandits et apprentissage par renforcement
# Implemente les equations de derivations/43_bandits_rl.qmd. R base.
# Apprendre en AGISSANT : arbitrer EXPLORATION (essayer pour apprendre) et
# EXPLOITATION (choisir le meilleur connu). Bandits (etat unique) puis MDP
# (etats, transitions) resolu par Q-learning.
# =============================================================================

#' Bandit UCB1 (borne de confiance superieure)
#'
#' Choisit a chaque tour le bras maximisant \eqn{\hat\mu_a+\sqrt{2\log t/n_a}} :
#' l'optimisme face a l'incertitude. Regret **logarithmique** \eqn{O(\log T)}.
#'
#' @param means moyennes VRAIES des bras (Bernoulli ou bornees)
#' @param horizon T
#' @param seed graine.
#' @return liste : `regret` (cumule), `arms` (choisis), `counts`.
#' @export
bandit_ucb <- function(means, horizon, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- length(means); n <- rep(0, K); s <- rep(0, K); arms <- integer(horizon)
  regret <- numeric(horizon); best <- max(means); cum <- 0
  for (t in seq_len(horizon)) {
    if (t <= K) a <- t else a <- which.max(s / n + sqrt(2 * log(t) / n))
    r <- rbinom(1, 1, means[a]); n[a] <- n[a] + 1; s[a] <- s[a] + r
    cum <- cum + (best - means[a]); regret[t] <- cum; arms[t] <- a
  }
  list(regret = regret, arms = arms, counts = n)
}

#' Bandit par echantillonnage de Thompson (Bernoulli)
#'
#' Prior Beta(1,1) par bras ; a chaque tour, tire \eqn{\theta_a\sim\text{Beta}
#' (\alpha_a,\beta_a)} et joue \eqn{\arg\max\theta_a} ; met a jour la posterieure.
#' Probability matching bayesien, regret logarithmique.
#'
#' @param means,horizon,seed cf. `bandit_ucb`.
#' @return liste : `regret`, `arms`, `counts`.
#' @export
bandit_thompson <- function(means, horizon, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- length(means); al <- rep(1, K); be <- rep(1, K); arms <- integer(horizon)
  regret <- numeric(horizon); best <- max(means); cum <- 0
  for (t in seq_len(horizon)) {
    a <- which.max(rbeta(K, al, be))
    r <- rbinom(1, 1, means[a]); al[a] <- al[a] + r; be[a] <- be[a] + (1 - r)
    cum <- cum + (best - means[a]); regret[t] <- cum; arms[t] <- a
  }
  list(regret = regret, arms = arms, counts = al + be - 2)
}

#' Iteration sur les valeurs (MDP a modele connu) — la reference optimale
#'
#' Itere l'operateur de Bellman \eqn{Q(s,a)=R(s,a)+\gamma\sum_{s'}P(s'|s,a)
#' \max_{a'}Q(s',a')} jusqu'a convergence. Fournit le \eqn{Q^\*} optimal.
#'
#' @param P transitions, tableau (S x A x S)
#' @param R recompenses (S x A)
#' @param gamma facteur d'actualisation
#' @param tol tolerance.
#' @return liste : `Q`, `V`, `policy`.
#' @export
value_iteration <- function(P, R, gamma = 0.9, tol = 1e-10) {
  S <- dim(P)[1]; A <- dim(P)[2]; Q <- matrix(0, S, A)
  repeat {
    V <- apply(Q, 1, max)
    Qn <- R + gamma * apply(P, c(1, 2), function(p) sum(p * V))
    if (max(abs(Qn - Q)) < tol) { Q <- Qn; break }; Q <- Qn
  }
  list(Q = Q, V = apply(Q, 1, max), policy = apply(Q, 1, which.max))
}

#' Q-learning tabulaire (MDP a modele INCONNU)
#'
#' Apprend \eqn{Q^\*} par interaction, sans connaitre \eqn{P,R} : politique
#' \eqn{\varepsilon}-greedy, mise a jour par difference temporelle
#' \eqn{Q(s,a)\leftarrow Q(s,a)+\alpha[r+\gamma\max_{a'}Q(s',a')-Q(s,a)]}.
#' Converge vers l'optimum du `value_iteration`.
#'
#' @param P,R,gamma cf. `value_iteration` (utilises pour SIMULER l'environnement) ;
#' @param episodes,steps duree
#' @param alpha,epsilon taux et exploration
#' @param seed graine.
#' @return liste : `Q`, `policy`.
#' @export
q_learning <- function(P, R, gamma = 0.9, episodes = 3000L, steps = 30L,
                       alpha = 0.1, epsilon = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  S <- dim(P)[1]; A <- dim(P)[2]; Q <- matrix(0, S, A)
  for (ep in seq_len(episodes)) {
    s <- sample.int(S, 1)
    for (t in seq_len(steps)) {
      a <- if (runif(1) < epsilon) sample.int(A, 1) else which.max(Q[s, ])
      s2 <- sample.int(S, 1, prob = P[s, a, ]); r <- R[s, a]
      Q[s, a] <- Q[s, a] + alpha * (r + gamma * max(Q[s2, ]) - Q[s, a])
      s <- s2
    }
  }
  list(Q = Q, policy = apply(Q, 1, which.max))
}
