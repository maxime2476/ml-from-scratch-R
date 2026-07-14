# =============================================================================
# Module 38 — Attention et Transformer minimal
# Implemente les equations de derivations/38_attention.qmd. R base.
# L'attention remplace la recurrence (Module 37) par un acces DIRECT a toutes les
# positions : chaque sortie est une moyenne ponderee des valeurs, les poids etant
# la SIMILARITE requete-cle. Base des Transformers (LLM, vision moderne).
# =============================================================================

#' Softmax numeriquement stable (par ligne)
#'
#' @param X matrice ; @param axis 1 (par ligne, defaut).
#' @return matrice de memes dimensions, chaque ligne sommant a 1.
#' @export
softmax_rows <- function(X) {
  X <- X - apply(X, 1, max); E <- exp(X); E / rowSums(E)
}

#' Attention par produit scalaire mis a l'echelle
#'
#' \eqn{\mathrm{Attn}(Q,K,V)=\mathrm{softmax}\!\bigl(QK^\top/\sqrt{d_k}\bigr)V}.
#' Chaque requete \eqn{q_i} interroge toutes les cles ; les poids (similarites
#' normalisees) melangent les valeurs. Une option de masque **causal** interdit
#' de regarder le futur (auto-regression).
#'
#' @param Q,K,V matrices requetes (T_q x d_k), cles (T_k x d_k), valeurs (T_k x d_v).
#' @param mask logique : masque causal (position i ne voit que j <= i).
#' @return liste : `out` (T_q x d_v), `weights` (T_q x T_k).
#' @export
attention <- function(Q, K, V, mask = FALSE) {
  dk <- ncol(K); scores <- (Q %*% t(K)) / sqrt(dk)
  if (mask) { Tq <- nrow(Q); m <- upper.tri(matrix(0, Tq, ncol(scores))); scores[m] <- -Inf }
  W <- softmax_rows(scores)
  list(out = W %*% V, weights = W)
}

#' Attention multi-tetes
#'
#' Projette \eqn{X} en \eqn{Q,K,V}, decoupe en `n_heads` sous-espaces, applique
#' l'attention en parallele dans chacun, concatene, puis reprojette. Chaque tete
#' capte un type de relation different.
#'
#' @param X entree (T x d_model).
#' @param Wq,Wk,Wv,Wo matrices de projection (d_model x d_model).
#' @param n_heads nombre de tetes (divise d_model) ; @param mask masque causal.
#' @return liste : `out` (T x d_model), `weights` (liste par tete).
#' @export
multi_head_attention <- function(X, Wq, Wk, Wv, Wo, n_heads = 1L, mask = FALSE) {
  d <- ncol(X); dh <- d %/% n_heads
  Q <- X %*% Wq; K <- X %*% Wk; V <- X %*% Wv
  heads <- vector("list", n_heads); wts <- vector("list", n_heads)
  for (h in seq_len(n_heads)) {
    idx <- ((h - 1) * dh + 1):(h * dh)
    a <- attention(Q[, idx, drop = FALSE], K[, idx, drop = FALSE], V[, idx, drop = FALSE], mask)
    heads[[h]] <- a$out; wts[[h]] <- a$weights
  }
  concat <- do.call(cbind, heads)
  list(out = concat %*% Wo, weights = wts)
}

#' Encodage positionnel sinusoidal
#'
#' \eqn{PE_{pos,2i}=\sin(pos/10000^{2i/d})}, \eqn{PE_{pos,2i+1}=\cos(\cdot)}.
#' Ajoute l'information d'ORDRE (que l'attention, permutation-equivariante, ignore).
#'
#' @param seq_len longueur de la sequence ; @param d_model dimension du modele.
#' @return matrice (seq_len x d_model).
#' @export
positional_encoding <- function(seq_len, d_model) {
  PE <- matrix(0, seq_len, d_model); pos <- 0:(seq_len - 1)
  for (i in seq_len(d_model %/% 2)) {
    div <- 10000^((2 * (i - 1)) / d_model)
    PE[, 2 * i - 1] <- sin(pos / div); PE[, 2 * i] <- cos(pos / div)
  }
  PE
}

#' Normalisation de couche (layer norm)
#'
#' Normalise chaque LIGNE (par exemple) puis remet a l'echelle : brique des blocs
#' de Transformer (avec les connexions residuelles).
#'
#' @param X matrice ; @param gamma,beta echelle et decalage (longueur ncol) ;
#' @param eps stabilisateur.
#' @return matrice normalisee.
#' @export
layer_norm <- function(X, gamma = rep(1, ncol(X)), beta = rep(0, ncol(X)), eps = 1e-5) {
  mu <- rowMeans(X); xc <- X - mu; v <- rowMeans(xc^2)
  xhat <- xc / sqrt(v + eps)
  sweep(sweep(xhat, 2, gamma, "*"), 2, beta, "+")
}
