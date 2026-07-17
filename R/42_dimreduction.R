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
#' @param X matrice n x p
#' @param k nb de composantes
#' @param gamma echelle RBF.
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
#' @param X melange (n x d)
#' @param n_comp nb de composantes
#' @param iter,tol arret.
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
#' @param V matrice non negative n x m
#' @param k rang
#' @param iter iterations.
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

#' Voisinages conditionnels a perplexite fixee (eq. 42.4-42.5)
#'
#' Calibre un sigma_i PAR POINT pour que chaque ligne de P atteigne la
#' perplexite cible. L'entropie de P_i croit avec sigma_i (Prop. 42.1), donc la
#' cible est atteinte en un sigma_i unique, obtenu par dichotomie sur
#' beta_i = 1/(2 sigma_i^2).
#'
#' @param D2 matrice n x n des distances au carre.
#' @param perplexity perplexite cible (nombre effectif de voisins).
#' @param tol tolerance sur log2(perplexite).
#' @param max_iter nombre maximal de pas de dichotomie.
#' @return matrice n x n des \eqn{p_{j|i}} (lignes sommant a 1, diagonale nulle).
#' @keywords internal
.p_conditional_perplexity <- function(D2, perplexity, tol = 1e-5, max_iter = 50L) {
  n <- nrow(D2)
  logU <- log2(perplexity)                       # entropie cible H = log2(Perp)
  P <- matrix(0, n, n)

  for (i in seq_len(n)) {
    d <- D2[i, -i]
    beta <- 1                                    # beta = 1/(2 sigma^2)
    lo <- -Inf; hi <- Inf

    for (k in seq_len(max_iter)) {
      # p non normalise, stabilise : le max des -beta*d vaut 0
      w <- exp(-d * beta - max(-d * beta))
      s <- sum(w)
      if (s < 1e-300) { beta <- beta / 2; next } # sigma trop petit : tout s'annule
      p <- w / s
      H <- -sum(p * log2(pmax(p, 1e-300)))       # entropie en bits

      if (abs(H - logU) < tol) break
      if (H > logU) {                            # trop diffus -> reduire sigma
        lo <- beta
        beta <- if (is.infinite(hi)) beta * 2 else (beta + hi) / 2
      } else {                                   # trop concentre -> agrandir sigma
        hi <- beta
        beta <- if (is.infinite(lo)) beta / 2 else (beta + lo) / 2
      }
    }
    P[i, -i] <- p
  }
  P
}

#' t-SNE (visualisation preservant le voisinage) — version compacte
#'
#' Convertit les distances en probabilites de voisinage (gaussiennes en haute
#' dimension, Student-t en basse) et minimise la divergence KL par descente de
#' gradient. Preserve la structure LOCALE : les clusters se separent nettement.
#'
#' @param X matrice n x p
#' @param dims dimension de sortie (2)
#' @param perplexity nombre effectif de voisins vise (eq. 42.5) ; fixe un
#'   sigma_i PAR POINT par dichotomie. Ignore si `sigma` est fourni.
#' @param sigma largeur commune imposee a tous les points ; par defaut `NULL`,
#'   et les sigma_i sont calibres depuis `perplexity`.
#' @param iter,eta nb d'iterations, pas de la descente.
#' @return matrice n x dims (le plongement).
#' @export
tsne <- function(X, dims = 2L, perplexity = 30, iter = 500L, eta = 200, sigma = NULL) {
  X <- as.matrix(X); n <- nrow(X)
  D2 <- outer(rowSums(X^2), rowSums(X^2), "+") - 2 * X %*% t(X); D2 <- pmax(D2, 0)
  if (is.null(sigma)) {
    P <- .p_conditional_perplexity(D2, perplexity)          # eq. 42.4-42.5
  } else {
    P <- exp(-D2 / (2 * sigma^2)); diag(P) <- 0; P <- P / rowSums(P)
  }
  P <- (P + t(P)) / (2 * n); P <- pmax(P, 1e-12)            # symetrisation (42.6)
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
