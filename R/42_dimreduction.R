# =============================================================================
# Module 42 — Reduction de dimension non lineaire
# Implemente les equations de derivations/42_dimreduction.qmd. R base.
# La PCA (Module 11) est LINEAIRE. Trois extensions : la PCA a NOYAU (structure
# courbe), l'ICA (sources INDEPENDANTES, pas seulement decorrelees), la NMF
# (parts ADDITIVES non negatives), et t-SNE (visualisation preservant le voisinage).
# =============================================================================

#' Analyse en composantes principales a noyau (kernel PCA)
#'
#' PCA dans l'espace de caracteristiques d'un noyau (Module 27) : on centre la
#' matrice de noyau \eqn{K_c=HKH} (\eqn{H=I-\tfrac1n\mathbf 1\mathbf 1^\top}), on
#' la diagonalise, et l'on projette. Capte des directions **non lineaires**.
#'
#' @param X matrice n x p ; @param k nb de composantes ; @param gamma echelle RBF.
#' @return liste : `proj` (n x k), `lambda` (valeurs propres).
#' @export
kernel_pca <- function(X, k = 2L, gamma = 1) {
  X <- as.matrix(X); n <- nrow(X)
  D2 <- outer(rowSums(X^2), rowSums(X^2), "+") - 2 * X %*% t(X)
  K <- exp(-gamma * pmax(D2, 0)); H <- diag(n) - 1 / n; Kc <- H %*% K %*% H
  ev <- eigen(Kc, symmetric = TRUE); lam <- ev$values[seq_len(k)]
  a <- sweep(ev$vectors[, seq_len(k), drop = FALSE], 2, sqrt(pmax(lam, 1e-10)), "/")
  list(proj = Kc %*% a, lambda = lam)
}

#' Analyse en composantes independantes (FastICA)
#'
#' Separe un melange \eqn{X=SA^\top} en **sources independantes** \eqn{S}, en
#' maximisant la **non-gaussianite** (negentropie, contraste \eqn{g=\tanh}) par
#' iteration de point fixe, apres **blanchiment**. Va au-dela de la PCA (qui ne
#' fait que decorreler) : resout le "cocktail party".
#'
#' @param X melange (n x d) ; @param n_comp nb de composantes ; @param iter,tol arret.
#' @return liste : `S` (sources estimees), `W` (matrice de separation).
#' @export
ica_fastica <- function(X, n_comp = ncol(X), iter = 300L, tol = 1e-9) {
  X <- scale(as.matrix(X), scale = FALSE); e <- eigen(cov(X))
  Wh <- e$vectors %*% diag(1 / sqrt(pmax(e$values, 1e-12))) %*% t(e$vectors)
  Z <- X %*% Wh                                             # blanchiment
  W <- matrix(rnorm(n_comp^2), n_comp, n_comp); W <- W / sqrt(rowSums(W^2))
  for (it in seq_len(iter)) {
    Wn <- matrix(0, n_comp, n_comp)
    for (p in seq_len(n_comp)) {
      wp <- W[p, ]; g <- tanh(Z %*% wp); gp <- 1 - g^2
      wp <- colMeans(Z * as.numeric(g)) - mean(gp) * wp
      if (p > 1) for (j in seq_len(p - 1)) wp <- wp - sum(wp * Wn[j, ]) * Wn[j, ]
      Wn[p, ] <- wp / sqrt(sum(wp^2))
    }
    if (max(abs(abs(rowSums(Wn * W)) - 1)) < tol) { W <- Wn; break }
    W <- Wn
  }
  list(S = Z %*% t(W), W = W %*% t(Wh))
}

#' Factorisation en matrices non negatives (NMF)
#'
#' \eqn{V\approx WH} avec \eqn{W,H\ge 0} (mises a jour multiplicatives de
#' Lee-Seung). La contrainte de non-negativite donne une decomposition en **parts
#' additives** interpretables (contrairement a la PCA, qui peut soustraire).
#'
#' @param V matrice non negative n x m ; @param k rang ; @param iter iterations.
#' @return liste : `W` (n x k), `H` (k x m), `reconstruction`.
#' @export
nmf <- function(V, k, iter = 500L) {
  V <- as.matrix(V); W <- matrix(runif(nrow(V) * k), nrow(V), k); H <- matrix(runif(k * ncol(V)), k, ncol(V))
  for (it in seq_len(iter)) {
    H <- H * (t(W) %*% V) / ((t(W) %*% W %*% H) + 1e-10)
    W <- W * (V %*% t(H)) / ((W %*% H %*% t(H)) + 1e-10)
  }
  list(W = W, H = H, reconstruction = W %*% H)
}

#' t-SNE (visualisation preservant le voisinage) — version compacte
#'
#' Convertit les distances en probabilites de voisinage (gaussiennes en haute
#' dimension, Student-t en basse) et minimise la divergence KL par descente de
#' gradient. Preserve la structure LOCALE : les clusters se separent nettement.
#'
#' @param X matrice n x p ; @param dims dimension de sortie (2) ; @param perplexity ;
#' @param iter,eta,sigma nb d'iterations, pas, largeur des gaussiennes.
#' @return matrice n x dims (le plongement).
#' @export
tsne <- function(X, dims = 2L, perplexity = 30, iter = 500L, eta = 200, sigma = NULL) {
  X <- as.matrix(X); n <- nrow(X)
  D2 <- outer(rowSums(X^2), rowSums(X^2), "+") - 2 * X %*% t(X); D2 <- pmax(D2, 0)
  if (is.null(sigma)) sigma <- sqrt(median(D2[D2 > 0]))
  P <- exp(-D2 / (2 * sigma^2)); diag(P) <- 0; P <- P / rowSums(P)
  P <- (P + t(P)) / (2 * n); P <- pmax(P, 1e-12)            # symetrisation
  set.seed(1); Y <- matrix(rnorm(n * dims) * 1e-4, n, dims); iY <- matrix(0, n, dims)
  for (t in seq_len(iter)) {
    dY <- outer(rowSums(Y^2), rowSums(Y^2), "+") - 2 * Y %*% t(Y)
    num <- 1 / (1 + pmax(dY, 0)); diag(num) <- 0; Q <- pmax(num / sum(num), 1e-12)
    L <- (P - Q) * num
    grad <- 4 * (diag(rowSums(L)) - L) %*% Y
    iY <- 0.8 * iY - eta * grad; Y <- Y + iY
    Y <- sweep(Y, 2, colMeans(Y))
  }
  Y
}
