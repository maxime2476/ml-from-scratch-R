# =============================================================================
# Module 36 — Reseaux de neurones convolutifs (CNN)
# Implemente les equations de derivations/36_cnn.qmd. R base.
# La convolution remplace la connexion dense par un PARTAGE DE POIDS local : un
# petit noyau balaie l'image, detectant le meme motif partout (equivariance a la
# translation). C'est la brique des reseaux de vision. Entree mono-canal ici.
# =============================================================================

#' Convolution 2D (correlation croisee, mode "valid")
#'
#' Pour chaque filtre \eqn{f}, \eqn{Y_{ij}^f=b_f+\sum_{a,b}X_{i+a-1,j+b-1}K_{ab}^f}.
#' Partage de poids : le meme noyau \eqn{K^f} est applique en toute position.
#'
#' @param X image d'entree (matrice H x W).
#' @param K noyaux, tableau (kh x kw x F)
#' @param b biais (longueur F).
#' @return liste : `out` (Hout x Wout x F), `cache`.
#' @export
conv2d <- function(X, K, b = rep(0, dim(K)[3])) {
  H <- nrow(X); W <- ncol(X); kh <- dim(K)[1]; kw <- dim(K)[2]; F <- dim(K)[3]
  oh <- H - kh + 1; ow <- W - kw + 1; out <- array(0, c(oh, ow, F))
  for (f in seq_len(F)) for (i in seq_len(oh)) for (j in seq_len(ow))
    out[i, j, f] <- b[f] + sum(X[i:(i + kh - 1), j:(j + kw - 1)] * K[, , f])
  list(out = out, cache = list(X = X, K = K, kh = kh, kw = kw, F = F, oh = oh, ow = ow))
}

#' Retropropagation de la convolution 2D
#'
#' Gradients par rapport a l'entree, aux noyaux et aux biais, etant donne `dout`.
#'
#' @param dout gradient en sortie (Hout x Wout x F).
#' @param cache sortie de `conv2d`.
#' @return liste : `dX`, `dK`, `db`.
#' @export
conv2d_backward <- function(dout, cache) {
  X <- cache$X; K <- cache$K; kh <- cache$kh; kw <- cache$kw; F <- cache$F
  dX <- matrix(0, nrow(X), ncol(X)); dK <- array(0, dim(K)); db <- numeric(F)
  for (f in seq_len(F)) for (i in seq_len(cache$oh)) for (j in seq_len(cache$ow)) {
    g <- dout[i, j, f]; db[f] <- db[f] + g
    patch <- i:(i + kh - 1); patchj <- j:(j + kw - 1)
    dK[, , f] <- dK[, , f] + g * X[patch, patchj]
    dX[patch, patchj] <- dX[patch, patchj] + g * K[, , f]
  }
  list(dX = dX, dK = dK, db = db)
}

#' Max-pooling 2D (sous-echantillonnage par le maximum)
#'
#' Reduit chaque bloc `pool` x `pool` a son maximum (invariance locale a la
#' translation, reduction de dimension).
#'
#' @param X carte de caracteristiques (H x W)
#' @param pool taille du bloc.
#' @return liste : `out`, `cache` (positions des maxima).
#' @export
max_pool2d <- function(X, pool = 2L) {
  H <- nrow(X); W <- ncol(X); oh <- H %/% pool; ow <- W %/% pool
  out <- matrix(0, oh, ow); argmax <- matrix(0, oh * ow, 2)
  k <- 0
  for (i in seq_len(oh)) for (j in seq_len(ow)) {
    ri <- ((i - 1) * pool + 1):((i - 1) * pool + pool); rj <- ((j - 1) * pool + 1):((j - 1) * pool + pool)
    blk <- X[ri, rj]; wm <- which.max(blk); out[i, j] <- blk[wm]
    k <- k + 1; argmax[k, ] <- c(ri[((wm - 1) %% pool) + 1], rj[((wm - 1) %/% pool) + 1])
  }
  list(out = out, cache = list(argmax = argmax, dimX = c(H, W), oh = oh, ow = ow))
}

#' Retropropagation du max-pooling
#'
#' Route le gradient vers la position du maximum de chaque bloc (les autres
#' recoivent 0).
#'
#' @param dout gradient en sortie (oh x ow)
#' @param cache sortie de `max_pool2d`.
#' @return `dX` (matrice de la taille de l'entree).
#' @export
max_pool2d_backward <- function(dout, cache) {
  dX <- matrix(0, cache$dimX[1], cache$dimX[2]); am <- cache$argmax; k <- 0
  for (i in seq_len(cache$oh)) for (j in seq_len(cache$ow)) {
    k <- k + 1; dX[am[k, 1], am[k, 2]] <- dX[am[k, 1], am[k, 2]] + dout[i, j]
  }
  dX
}
