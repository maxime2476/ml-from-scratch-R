# =============================================================================
# Module 39 — Machines a vecteurs de support (SVM)
# Implemente les equations de derivations/39_svm.qmd. R base + quadprog.
# On cherche l'hyperplan de MARGE MAXIMALE separant deux classes. La formulation
# DUALE ne fait intervenir les donnees que par produits scalaires -> l'astuce du
# NOYAU permet des frontieres non lineaires sans jamais expliciter l'espace.
# =============================================================================

#' Noyau gaussien (RBF) pour SVM
#' @param X1,X2 matrices
#' @param gamma largeur inverse.
#' @return matrice de noyau.
#' @export
svm_rbf <- function(X1, X2, gamma) {
  D2 <- outer(rowSums(X1^2), rowSums(X2^2), "+") - 2 * X1 %*% t(X2)
  exp(-gamma * pmax(D2, 0))
}

#' SVM a marge souple par le probleme dual (programmation quadratique)
#'
#' Maximise \eqn{\sum_i\alpha_i-\tfrac12\sum_{ij}\alpha_i\alpha_j y_iy_j K(x_i,x_j)}
#' sous \eqn{0\le\alpha_i\le C} et \eqn{\sum_i\alpha_i y_i=0}. Resolu par
#' `quadprog::solve.QP`. Seuls les **vecteurs de support** (\eqn{\alpha_i>0})
#' comptent.
#'
#' @param X matrice n x p
#' @param y etiquettes dans \eqn{\{-1,+1\}}.
#' @param C penalite (marge souple)
#' @param kernel "linear" ou "rbf"
#' @param gamma parametre RBF.
#' @return objet `svm` : `alpha`, `b`, vecteurs de support, hyperparametres.
#' @export
svm_fit <- function(X, y, C = 1, kernel = c("linear", "rbf"), gamma = 0.5) {
  kernel <- match.arg(kernel); X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X)
  Kmat <- if (kernel == "linear") X %*% t(X) else svm_rbf(X, X, gamma)
  Dmat <- (y %o% y) * Kmat + diag(1e-8, n)                 # rendre definie positive
  Amat <- cbind(y, diag(n), -diag(n)); bvec <- c(0, rep(0, n), rep(-C, n))
  sol <- quadprog::solve.QP(Dmat, rep(1, n), Amat, bvec, meq = 1)
  a <- sol$solution; a[a < 1e-6] <- 0
  margin_sv <- which(a > 1e-6 & a < C - 1e-6)
  f_nob <- as.numeric((a * y) %*% Kmat)
  b <- if (length(margin_sv)) mean(y[margin_sv] - f_nob[margin_sv]) else 0
  structure(list(alpha = a, b = b, y = y, X = X, kernel = kernel, gamma = gamma,
                 C = C, sv = which(a > 1e-6), n_sv = sum(a > 1e-6)), class = "svm_scratch")
}

#' Prediction d'un SVM
#'
#' \eqn{\hat y(x)=\mathrm{sign}\bigl(\sum_i\alpha_i y_i K(x_i,x)+b\bigr)}.
#'
#' @param object objet `svm_fit`
#' @param Xnew nouvelles observations
#' @param decision renvoyer la valeur de decision (defaut FALSE : le signe).
#' @return vecteur de classes (ou de valeurs de decision).
#' @export
svm_predict <- function(object, Xnew, decision = FALSE) {
  Xnew <- as.matrix(Xnew)
  Kx <- if (object$kernel == "linear") Xnew %*% t(object$X) else svm_rbf(Xnew, object$X, object$gamma)
  f <- as.numeric(Kx %*% (object$alpha * object$y)) + object$b
  if (decision) f else sign(f)
}
