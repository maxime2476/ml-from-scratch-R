# =============================================================================
# Module 11 — Non supervisé : ACP, k-means, EM (mélange gaussien)
# Implémente les équations de derivations/11_pca_kmeans_em.qmd. R base + Module 0.
# =============================================================================

#' Analyse en composantes principales par SVD (éq. 11.2-11.3)
#'
#' Centre X puis calcule l'ACP via la SVD (voie 2 = voie 1 par Prop. 11.2).
#' Reproduit `prcomp` : `sdev` = d/sqrt(n-1), `rotation` = vecteurs singuliers
#' droits V, `scores` = X_centré V.
#'
#' @param X matrice n x p.
#' @param center centrer les colonnes (défaut TRUE).
#' @param scale réduire les colonnes (défaut FALSE).
#' @return liste : `sdev`, `rotation` (p x p), `scores` (n x p), `var_explained`,
#'   `center`, `scale`.
#' @export
pca_fit <- function(X, center = TRUE, scale = FALSE) {
  X <- as.matrix(X); n <- nrow(X)
  ctr <- if (center) colMeans(X) else rep(0, ncol(X))
  Xc <- sweep(X, 2, ctr, "-")
  scl <- if (scale) apply(Xc, 2, sd) else rep(1, ncol(X))
  Xs <- sweep(Xc, 2, scl, "/")
  sv <- svd(Xs)
  sdev <- sv$d / sqrt(n - 1)
  scores <- Xs %*% sv$v
  list(sdev = sdev, rotation = sv$v, scores = scores,
       var_explained = sdev^2 / sum(sdev^2), center = ctr, scale = scl)
}

#' k-means par l'algorithme de Lloyd (éq. 11.4)
#'
#' Minimisation alternée de l'inertie intra-classe (affectation au centre le plus
#' proche, puis moyenne de classe). Reproduit `stats::kmeans(algorithm="Lloyd")`
#' à initialisation identique.
#'
#' @param X matrice n x p.
#' @param K nombre de classes.
#' @param centers matrice K x p de centres initiaux (sinon tirage aléatoire).
#' @param max_iter itérations maximales.
#' @param nstart nombre de redémarrages aléatoires (si centers non fourni).
#' @param seed graine.
#' @return liste : `cluster`, `centers`, `withinss`, `tot_withinss`, `iter`.
#' @export
kmeans_fit <- function(X, K, centers = NULL, max_iter = 100L, nstart = 1L, seed = NULL) {
  X <- as.matrix(X); n <- nrow(X)
  if (!is.null(seed)) set.seed(seed)
  run_once <- function(cen) {
    cl <- integer(n); it <- 0L
    for (it in seq_len(max_iter)) {
      D <- sapply(seq_len(K), function(k) colSums((t(X) - cen[k, ])^2))  # n x K
      new_cl <- max.col(-D, ties.method = "first")                       # plus proche centre
      cen <- t(sapply(seq_len(K), function(k)
        if (any(new_cl == k)) colMeans(X[new_cl == k, , drop = FALSE]) else cen[k, ]))
      if (identical(new_cl, cl)) { cl <- new_cl; break }
      cl <- new_cl
    }
    wss <- sapply(seq_len(K), function(k) sum((t(X[cl == k, , drop = FALSE]) - cen[k, ])^2))
    list(cluster = cl, centers = cen, withinss = wss, tot_withinss = sum(wss), iter = it)
  }
  starts <- if (is.null(centers)) lapply(seq_len(nstart), function(s) X[sample.int(n, K), , drop = FALSE])
            else list(as.matrix(centers))
  best <- NULL
  for (cen in starts) {
    res <- run_once(cen)
    if (is.null(best) || res$tot_withinss < best$tot_withinss) best <- res
  }
  best
}

# Densité gaussienne multivariée (log), stable par Cholesky.
.log_dmvnorm <- function(X, mu, Sigma) {
  p <- length(mu)
  R <- chol(Sigma)                                   # Sigma = R'R (R triangulaire sup.)
  Xc <- sweep(X, 2, mu, "-")
  z <- backsolve(R, t(Xc), transpose = TRUE)          # z = R^{-T}(x-mu)
  quad <- colSums(z^2)
  logdet <- 2 * sum(log(diag(R)))
  -0.5 * (p * log(2 * pi) + logdet + quad)
}

# log-sum-exp par ligne (stabilité numérique).
.logsumexp_rows <- function(M) {
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

#' Log-vraisemblance observée d'un mélange gaussien
#'
#' @param X données n x p.
#' @param pi poids (longueur K).
#' @param mu liste (ou matrice K x p) des moyennes.
#' @param Sigma liste des K covariances.
#' @return la log-vraisemblance observée.
#' @export
gmm_loglik <- function(X, pi, mu, Sigma) {
  X <- as.matrix(X); K <- length(pi)
  if (is.matrix(mu)) mu <- lapply(seq_len(K), function(k) mu[k, ])
  logd <- sapply(seq_len(K), function(k) log(pi[k]) + .log_dmvnorm(X, mu[[k]], Sigma[[k]]))
  sum(.logsumexp_rows(logd))
}

#' Mélange gaussien par EM (éq. 11.5-11.7)
#'
#' Étape E : responsabilités (éq. 11.6) ; étape M : mises à jour fermées (éq.
#' 11.7). Initialisation par k-means. Convergence sur la log-vraisemblance.
#'
#' @param X matrice n x p.
#' @param K nombre de composantes.
#' @param max_iter itérations maximales.
#' @param tol tolérance sur la variation de log-vraisemblance.
#' @param reg régularisation ajoutée à la diagonale des covariances.
#' @param seed graine (initialisation k-means).
#' @return liste : `pi`, `mu` (K x p), `Sigma` (liste), `gamma` (responsabilités),
#'   `loglik`, `iter`, `cluster`.
#' @export
em_gmm <- function(X, K, max_iter = 200L, tol = 1e-8, reg = 1e-6, seed = NULL) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  km <- kmeans_fit(X, K, seed = seed, nstart = 5L)
  cl <- km$cluster
  pi <- as.numeric(table(factor(cl, levels = seq_len(K)))) / n
  mu <- km$centers
  Sigma <- lapply(seq_len(K), function(k) {
    Xk <- X[cl == k, , drop = FALSE]
    S <- if (nrow(Xk) > 1) cov(Xk) else diag(p)
    S + reg * diag(p)
  })

  loglik_old <- -Inf
  for (iter in seq_len(max_iter)) {
    # --- Étape E (éq. 11.6) ---
    logd <- sapply(seq_len(K), function(k) log(pi[k]) + .log_dmvnorm(X, mu[k, ], Sigma[[k]]))
    ls <- .logsumexp_rows(logd)
    loglik <- sum(ls)
    gamma <- exp(logd - ls)                          # responsabilités n x K
    # --- Étape M (éq. 11.7) ---
    Nk <- colSums(gamma)
    pi <- Nk / n
    mu <- (t(gamma) %*% X) / Nk
    Sigma <- lapply(seq_len(K), function(k) {
      Xc <- sweep(X, 2, mu[k, ], "-")
      (t(Xc * gamma[, k]) %*% Xc) / Nk[k] + reg * diag(p)
    })
    if (abs(loglik - loglik_old) < tol) break
    loglik_old <- loglik
  }
  list(pi = pi, mu = mu, Sigma = Sigma, gamma = gamma, loglik = loglik,
       iter = iter, cluster = max.col(gamma, ties.method = "first"))
}
