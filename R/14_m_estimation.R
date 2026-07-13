# =============================================================================
# Module 14 — M-estimation et lecture bayésienne (helpers)
# Implémente les objets de derivations/14_m_estimation.qmd. Module théorique :
# ces fonctions servent aux vérifications et illustrations.
# =============================================================================

#' Variance sandwich d'un M-estimateur (éq. 14.2)
#'
#' \eqn{\hat V = \hat A^{-1}\hat B\hat A^{-1}/n}, avec \eqn{\hat B = \frac1n\sum_i
#' \psi_i\psi_i^\top} et \eqn{\hat A} fournie (moyenne des dérivées de la fonction
#' d'estimation, \eqn{A = \mathbb E[-\partial_\theta\psi]}).
#'
#' @param psi matrice n x p des valeurs de la fonction d'estimation \eqn{\psi_i}.
#' @param A matrice p x p, \eqn{\hat A} (p.ex. \eqn{X^\top X/n} pour l'OLS).
#' @return la matrice de variance sandwich p x p de \eqn{\hat\theta}.
#' @export
m_estimation_vcov <- function(psi, A) {
  psi <- as.matrix(psi); n <- nrow(psi)
  B <- crossprod(psi) / n                 # (1/n) sum psi psi'
  Ainv <- solve(A)
  Ainv %*% B %*% t(Ainv) / n
}

#' Fonction d'estimation de l'OLS : psi_i = x_i * e_i
#'
#' @param X matrice de design n x p (constante incluse).
#' @param resid résidus \eqn{\hat\varepsilon_i}.
#' @return matrice n x p des \eqn{\psi_i = x_i \hat\varepsilon_i}.
#' @export
ols_psi <- function(X, resid) as.matrix(X) * resid

#' Postérieure conjuguée du ridge (Prop. 14.2)
#'
#' Sous \eqn{y\mid\beta\sim\mathcal N(X\beta,\sigma^2 I)} et prior
#' \eqn{\beta\sim\mathcal N(0,\tau^2 I)}, renvoie la moyenne et la covariance a
#' posteriori (éq. 14.3), avec \eqn{\lambda=\sigma^2/\tau^2}.
#'
#' @param X matrice de design n x p (sans intercept, ou déjà centré).
#' @param y réponse.
#' @param lambda pénalité ridge \eqn{\lambda=\sigma^2/\tau^2}.
#' @param sigma2 variance du bruit (connue).
#' @return liste : `mean` (= estimateur ridge), `cov` (\eqn{\Sigma_{\text{post}}}).
#' @export
ridge_posterior <- function(X, y, lambda, sigma2) {
  X <- as.matrix(X); y <- as.numeric(y); p <- ncol(X)
  Sigma_post <- sigma2 * solve(crossprod(X) + lambda * diag(p))
  mu_post <- as.numeric(solve(crossprod(X) + lambda * diag(p), crossprod(X, y)))
  list(mean = mu_post, cov = Sigma_post)
}
