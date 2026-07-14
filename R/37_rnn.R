# =============================================================================
# Module 37 — Reseaux recurrents (RNN, LSTM)
# Implemente les equations de derivations/37_rnn.qmd. R base.
# Pour les SEQUENCES, un etat cache h_t resume le passe et se propage dans le
# temps. L'apprentissage se fait par retropropagation DANS LE TEMPS (BPTT). Le
# RNN simple souffre du gradient qui s'evanouit ; la LSTM y remedie par une
# "autoroute" de memoire (l'etat de cellule).
# =============================================================================

.sigmoid <- function(x) 1 / (1 + exp(-x))

#' Passe avant d'un RNN simple (sequence -> sequence)
#'
#' \eqn{h_t=\tanh(W_{xh}x_t+W_{hh}h_{t-1}+b_h)}, \eqn{y_t=W_{hy}h_t+b_y}.
#'
#' @param X sequence d'entree (T x d).
#' @param Wxh,Whh,Why matrices de poids (H x d), (H x H), (O x H).
#' @param bh,by biais (H), (O).
#' @param h0 etat initial (defaut 0).
#' @return liste : `Y` (T x O), `H` (T x H), `cache`.
#' @export
rnn_forward <- function(X, Wxh, Whh, Why, bh, by, h0 = rep(0, nrow(Wxh))) {
  Tn <- nrow(X); H <- nrow(Wxh); O <- nrow(Why)
  Hs <- matrix(0, Tn, H); Ys <- matrix(0, Tn, O); h <- h0
  for (t in seq_len(Tn)) {
    h <- tanh(as.numeric(Wxh %*% X[t, ] + Whh %*% h + bh)); Hs[t, ] <- h
    Ys[t, ] <- as.numeric(Why %*% h + by)
  }
  list(Y = Ys, H = Hs, cache = list(X = X, Hs = Hs, Wxh = Wxh, Whh = Whh, Why = Why, h0 = h0))
}

#' Retropropagation dans le temps (BPTT) du RNN simple
#'
#' @param dY gradient en sortie (T x O).
#' @param cache sortie de `rnn_forward`.
#' @return liste : `dWxh`, `dWhh`, `dWhy`, `dbh`, `dby`.
#' @export
rnn_backward <- function(dY, cache) {
  X <- cache$X; Hs <- cache$Hs; Whh <- cache$Whh; Why <- cache$Why; Tn <- nrow(X)
  dWxh <- matrix(0, nrow(cache$Wxh), ncol(cache$Wxh)); dWhh <- matrix(0, nrow(Whh), ncol(Whh))
  dWhy <- matrix(0, nrow(Why), ncol(Why)); dbh <- rep(0, nrow(Whh)); dby <- rep(0, nrow(Why))
  dh_next <- rep(0, nrow(Whh))
  for (t in rev(seq_len(Tn))) {
    ht <- Hs[t, ]; h_prev <- if (t > 1) Hs[t - 1, ] else cache$h0
    dWhy <- dWhy + dY[t, ] %o% ht; dby <- dby + dY[t, ]
    dh <- as.numeric(t(Why) %*% dY[t, ]) + dh_next
    da <- dh * (1 - ht^2)                                   # derivee de tanh
    dbh <- dbh + da; dWxh <- dWxh + da %o% X[t, ]; dWhh <- dWhh + da %o% h_prev
    dh_next <- as.numeric(t(Whh) %*% da)
  }
  list(dWxh = dWxh, dWhh = dWhh, dWhy = dWhy, dbh = dbh, dby = dby)
}

#' Passe avant d'une cellule LSTM (sequence)
#'
#' Portes d'entree/oubli/sortie et etat de cellule :
#' \eqn{i=\sigma(W_i[x_t,h_{t-1}])}, \eqn{f=\sigma(\cdot)}, \eqn{o=\sigma(\cdot)},
#' \eqn{g=\tanh(\cdot)} ; \eqn{c_t=f\odot c_{t-1}+i\odot g}, \eqn{h_t=o\odot\tanh(c_t)}.
#' L'etat de cellule \eqn{c_t} forme une "autoroute" ou le gradient circule sans
#' s'evanouir.
#'
#' @param X sequence (T x d).
#' @param Wi,Wf,Wo,Wg matrices (H x (d+H)) empilant \eqn{[x_t,h_{t-1}]}.
#' @param bi,bf,bo,bg biais (H).
#' @return liste : `H` (T x H), `C` (T x H).
#' @export
lstm_forward <- function(X, Wi, Wf, Wo, Wg, bi, bf, bo, bg) {
  Tn <- nrow(X); H <- nrow(Wi); h <- rep(0, H); c <- rep(0, H)
  Hs <- matrix(0, Tn, H); Cs <- matrix(0, Tn, H)
  for (t in seq_len(Tn)) {
    z <- c(X[t, ], h)
    i <- .sigmoid(as.numeric(Wi %*% z + bi)); f <- .sigmoid(as.numeric(Wf %*% z + bf))
    o <- .sigmoid(as.numeric(Wo %*% z + bo)); g <- tanh(as.numeric(Wg %*% z + bg))
    c <- f * c + i * g; h <- o * tanh(c); Hs[t, ] <- h; Cs[t, ] <- c
  }
  list(H = Hs, C = Cs)
}
