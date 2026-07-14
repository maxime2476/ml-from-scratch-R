# =============================================================================
# Module 27 — Processus gaussiens et méthodes à noyaux (RKHS)
# Implémente les équations de derivations/27_gaussian_process.qmd. R base.
# Pont bayésien <-> fréquentiste : la moyenne a posteriori d'un GP EST la
# régression ridge à noyau (théorème de représentation), et prolonge le ridge
# du Module 4 dans un espace de caractéristiques de dimension infinie.
# =============================================================================

#' Noyau gaussien (RBF)
#'
#' \eqn{k(x,x')=\sigma_f^2\exp(-\|x-x'\|^2/(2\ell^2))}.
#'
#' @param X1,X2 matrices (n1 x p), (n2 x p).
#' @param lengthscale échelle \eqn{\ell}.
#' @param variance variance du signal \eqn{\sigma_f^2}.
#' @return matrice de noyau n1 x n2.
#' @export
rbf_kernel <- function(X1, X2, lengthscale = 1, variance = 1) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2)
  D2 <- outer(rowSums(X1^2), rowSums(X2^2), "+") - 2 * X1 %*% t(X2)
  variance * exp(-pmax(D2, 0) / (2 * lengthscale^2))
}

#' Ajustement d'un processus gaussien (régression)
#'
#' Résout \eqn{(K+\sigma_n^2 I)\alpha=y} par Cholesky (Module 0). Renvoie de quoi
#' prédire moyenne et variance a posteriori, et la **vraisemblance marginale**.
#'
#' @param X,y données d'apprentissage.
#' @param lengthscale,sigma_f,sigma_n hyperparamètres (échelle, signal, bruit).
#' @return objet `gp` : `alpha`, `L`, `X`, hyperparamètres, `loglik`.
#' @export
gp_fit <- function(X, y, lengthscale = 1, sigma_f = 1, sigma_n = 0.1) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X)
  K <- rbf_kernel(X, X, lengthscale, sigma_f^2) + sigma_n^2 * diag(n)
  L <- chol(K)                                              # K = L'L (L triangulaire sup)
  alpha <- backsolve(L, backsolve(L, y, transpose = TRUE))
  # log-vraisemblance marginale (éq. 27.3)
  loglik <- -0.5 * sum(y * alpha) - sum(log(diag(L))) - 0.5 * n * log(2 * pi)
  structure(list(alpha = alpha, L = L, X = X, y = y,
                 lengthscale = lengthscale, sigma_f = sigma_f, sigma_n = sigma_n,
                 loglik = loglik), class = "gp")
}

#' Prédiction d'un processus gaussien
#'
#' Moyenne a posteriori \eqn{\bar f_*=k_*^\top\alpha} et variance de la **fonction
#' latente** \eqn{\operatorname{Var}[f_*]=k_{**}-k_*^\top(K+\sigma_n^2I)^{-1}k_*}.
#' La variance d'une **observation** bruitée ajoute \eqn{\sigma_n^2}.
#'
#' @param object objet `gp`.
#' @param Xnew points de prédiction.
#' @return liste : `mean`, `sd` (fonction latente), `sd_obs` (observation bruitée).
#' @export
gp_predict <- function(object, Xnew) {
  Xnew <- as.matrix(Xnew)
  Ks <- rbf_kernel(object$X, Xnew, object$lengthscale, object$sigma_f^2)
  kss <- rep(object$sigma_f^2, nrow(Xnew))                  # diag k(x*,x*)
  mean <- as.numeric(t(Ks) %*% object$alpha)
  v <- backsolve(object$L, Ks, transpose = TRUE)
  var <- pmax(kss - colSums(v^2), 0)
  list(mean = mean, sd = sqrt(var), sd_obs = sqrt(var + object$sigma_n^2))
}

#' Régression ridge à noyau (théorème de représentation)
#'
#' \eqn{\hat f(x)=k(x)^\top(K+\lambda I)^{-1}y}. **Identique** à la moyenne a
#' posteriori du GP avec \eqn{\lambda=\sigma_n^2} (pont bayésien/fréquentiste).
#'
#' @param X,y données ; @param lengthscale,sigma_f noyau ; @param lambda pénalité.
#' @return fonction `newX -> prédictions`.
#' @export
kernel_ridge <- function(X, y, lengthscale = 1, sigma_f = 1, lambda = 0.01) {
  X <- as.matrix(X); K <- rbf_kernel(X, X, lengthscale, sigma_f^2)
  alpha <- solve(K + lambda * diag(nrow(X)), y)
  function(newX) as.numeric(rbf_kernel(as.matrix(newX), X, lengthscale, sigma_f^2) %*% alpha)
}

#' Sélection des hyperparamètres par maximum de vraisemblance marginale
#'
#' Maximise la log-vraisemblance marginale (éq. 27.3) sur
#' \eqn{(\ell,\sigma_f,\sigma_n)} — le rasoir d'Occam bayésien automatique.
#' Optimisation sur l'échelle log (positivité) via `optim` (L-BFGS-B).
#'
#' @param X,y données ; @param init valeurs initiales (échelle, signal, bruit).
#' @return objet `gp` optimisé (avec les hyperparamètres retenus).
#' @export
gp_optimize <- function(X, y, init = c(1, 1, 0.1)) {
  nll <- function(lp) {
    p <- exp(lp)
    -gp_fit(X, y, lengthscale = p[1], sigma_f = p[2], sigma_n = p[3])$loglik
  }
  opt <- optim(log(init), nll, method = "L-BFGS-B")
  p <- exp(opt$par)
  gp_fit(X, y, lengthscale = p[1], sigma_f = p[2], sigma_n = p[3])
}
