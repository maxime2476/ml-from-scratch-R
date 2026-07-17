# =============================================================================
# Module 40 — Classifieurs generatifs (LDA, QDA, Naive Bayes)
# Implemente les equations de derivations/40_generative.qmd. R base.
# Au lieu de modeliser directement P(y|x) (discriminatif : GLM, SVM), on modelise
# la GENERATION P(x|y) et le prior P(y), puis on inverse par Bayes. Hypotheses
# gaussiennes -> frontieres lineaires (LDA), quadratiques (QDA) ou naives (NB).
# =============================================================================

.class_stats <- function(X, y) {
  y <- as.factor(y); ks <- levels(y)
  list(ks = ks, prior = as.numeric(table(y) / length(y)),
       mu = lapply(ks, function(k) colMeans(X[y == k, , drop = FALSE])), y = y)
}

#' Analyse discriminante lineaire (LDA)
#'
#' Modele : \eqn{x\mid y=k\sim\mathcal N(\mu_k,\Sigma)} (covariance **commune**).
#' La regle de Bayes donne une frontiere **lineaire** ; le discriminant est
#' \eqn{\delta_k(x)=x^\top\Sigma^{-1}\mu_k-\tfrac12\mu_k^\top\Sigma^{-1}\mu_k+\log\pi_k}.
#'
#' @param X matrice n x p
#' @param y etiquettes de classe.
#' @return objet `lda_scratch`.
#' @export
lda_fit <- function(X, y) {
  s <- .class_stats(X, y); X <- as.matrix(X); n <- nrow(X)
  Sw <- Reduce("+", lapply(seq_along(s$ks), function(k) {
    Xc <- sweep(X[s$y == s$ks[k], , drop = FALSE], 2, s$mu[[k]]); crossprod(Xc)
  })) / (n - length(s$ks))
  structure(c(s, list(Sinv = solve(Sw))), class = "lda_scratch")
}

#' Prediction LDA
#' @param object objet `lda_fit`
#' @param Xnew nouvelles observations.
#' @return classes predites.
#' @export
lda_predict <- function(object, Xnew) {
  Xnew <- as.matrix(Xnew)
  D <- sapply(seq_along(object$ks), function(k) { mu <- object$mu[[k]]
    as.numeric(Xnew %*% object$Sinv %*% mu) - 0.5 * as.numeric(t(mu) %*% object$Sinv %*% mu) + log(object$prior[k])
  })
  object$ks[max.col(D)]
}

#' Analyse discriminante quadratique (QDA)
#'
#' Chaque classe a sa **propre** covariance \eqn{\Sigma_k} : frontiere
#' **quadratique**. \eqn{\delta_k(x)=-\tfrac12\log|\Sigma_k|-\tfrac12(x-\mu_k)^\top
#' \Sigma_k^{-1}(x-\mu_k)+\log\pi_k}.
#'
#' @param X,y cf. `lda_fit`.
#' @return objet `qda_scratch`.
#' @export
qda_fit <- function(X, y) {
  s <- .class_stats(X, y); X <- as.matrix(X)
  Sig <- lapply(seq_along(s$ks), function(k) {
    Xc <- sweep(X[s$y == s$ks[k], , drop = FALSE], 2, s$mu[[k]]); crossprod(Xc) / (sum(s$y == s$ks[k]) - 1)
  })
  structure(c(s, list(Sigma = Sig)), class = "qda_scratch")
}

#' Prediction QDA
#' @param object objet `qda_fit`
#' @param Xnew nouvelles observations.
#' @return classes predites.
#' @export
qda_predict <- function(object, Xnew) {
  Xnew <- as.matrix(Xnew)
  D <- sapply(seq_along(object$ks), function(k) {
    mu <- object$mu[[k]]; Si <- solve(object$Sigma[[k]])
    ld <- as.numeric(determinant(object$Sigma[[k]], logarithm = TRUE)$modulus)
    Xc <- sweep(Xnew, 2, mu); -0.5 * ld - 0.5 * rowSums((Xc %*% Si) * Xc) + log(object$prior[k])
  })
  object$ks[max.col(D)]
}

#' Classifieur naif de Bayes (gaussien)
#'
#' Hypothese **naive** : les features sont independantes conditionnellement a la
#' classe, \eqn{P(x\mid y=k)=\prod_j\mathcal N(x_j;\mu_{kj},\sigma_{kj}^2)}. Rapide,
#' robuste en haute dimension, mais biaise si les features sont correlees.
#'
#' @param X,y cf. `lda_fit`.
#' @return objet `nb_scratch`.
#' @export
naive_bayes_fit <- function(X, y) {
  s <- .class_stats(X, y); X <- as.matrix(X)
  sds <- lapply(seq_along(s$ks), function(k) apply(X[s$y == s$ks[k], , drop = FALSE], 2, sd))
  structure(c(s, list(sd = sds)), class = "nb_scratch")
}

#' Prediction Naive Bayes
#' @param object objet `naive_bayes_fit`
#' @param Xnew nouvelles observations.
#' @return classes predites.
#' @export
naive_bayes_predict <- function(object, Xnew) {
  Xnew <- as.matrix(Xnew)
  D <- sapply(seq_along(object$ks), function(k) {
    ll <- rowSums(sapply(seq_len(ncol(Xnew)), function(j)
      dnorm(Xnew[, j], object$mu[[k]][j], object$sd[[k]][j], log = TRUE)))
    ll + log(object$prior[k])
  })
  object$ks[max.col(D)]
}
