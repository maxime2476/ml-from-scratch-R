# =============================================================================
# Module 24 — Fonctions d'influence : le fil unificateur
# Implémente les équations de derivations/24_influence.qmd. R base.
# Idée : TOUT estimateur asymptotiquement linéaire s'écrit
#   sqrt(n)(theta_hat - theta) = (1/sqrt(n)) sum_i IC_i + o_p(1),
# et sa variance, son bootstrap (M17), son jackknife et son sandwich (M2/M14)
# ne sont que des façons d'estimer Var(IC). Le score de Neyman (M16) et la
# correction du lasso débiaisé (M22) sont des fonctions d'influence EFFICACES.
# =============================================================================

#' Fonction d'influence de l'OLS
#'
#' \eqn{\mathrm{IC}_i = n\,(X^\top X)^{-1} x_i \hat\varepsilon_i}, de sorte que
#' \eqn{\hat\beta-\beta \approx \frac1n\sum_i \mathrm{IC}_i} et
#' \eqn{\widehat{\operatorname{Var}}(\hat\beta)=\frac1{n^2}\sum_i \mathrm{IC}_i\mathrm{IC}_i^\top}
#' — **exactement** le sandwich HC0 (Module 2).
#'
#' @param X matrice de design n x p (constante incluse).
#' @param y réponse.
#' @return liste : `ic` (n x p), `vcov` (= sandwich HC0), `beta`.
#' @export
influence_ols <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X)
  XtXinv <- solve(crossprod(X))
  beta <- as.numeric(XtXinv %*% crossprod(X, y))
  e <- as.numeric(y - X %*% beta)
  ic <- n * (X * e) %*% t(XtXinv)                 # ligne i = n (X'X)^{-1} x_i e_i
  list(ic = ic, vcov = crossprod(ic) / n^2, beta = beta)
}

#' Fonction d'influence d'un MLE (GLM canonique)
#'
#' \eqn{\mathrm{IC}_i = n\,\mathcal I^{-1} x_i (y_i-\hat\mu_i)}, avec
#' \eqn{\mathcal I = X^\top W X} l'information de Fisher (Module 3). Sa variance
#' \eqn{\frac1{n^2}\sum \mathrm{IC}_i\mathrm{IC}_i^\top} est le sandwich robuste ;
#' sous spécification correcte elle vaut \eqn{\mathcal I^{-1}} (variance modèle).
#'
#' @param fit objet `glm_irls` (Module 3, lien canonique).
#' @return liste : `ic` (n x p), `vcov`.
#' @export
influence_mle <- function(fit) {
  X <- fit$model_matrix; n <- nrow(X)
  s <- X * (fit$response - fit$fitted)            # score canonique x_i (y_i - mu_i)
  ic <- n * s %*% t(fit$vcov)                      # I^{-1} = fit$vcov
  list(ic = ic, vcov = crossprod(ic) / n^2)
}

#' Jackknife (delete-one) — le jackknife infinitésimal EST la fonction d'influence
#'
#' Pour un estimateur lisse, \eqn{\hat\theta_{(-i)}-\hat\theta \approx -\mathrm{IC}_i/n},
#' d'où \eqn{\widehat{\operatorname{Var}}_{\text{jack}} = \frac{n-1}{n}\sum_i
#' (\hat\theta_{(-i)}-\bar\theta_{(\cdot)})^2 \approx \frac1{n^2}\sum \mathrm{IC}_i^2}.
#'
#' @param data data.frame (une ligne = une observation).
#' @param stat_fn fonction `data -> theta` (scalaire).
#' @return liste : `estimate`, `bias`, `var`, `values`.
#' @export
jackknife <- function(data, stat_fn) {
  n <- nrow(data)
  th <- stat_fn(data)
  loo <- vapply(seq_len(n), function(i) stat_fn(data[-i, , drop = FALSE]), numeric(1))
  thbar <- mean(loo)
  list(estimate = th, bias = (n - 1) * (thbar - th),
       var = (n - 1) / n * sum((loo - thbar)^2), values = loo)
}

#' Estimateur « en un pas » (one-step / Newton-scoring)
#'
#' \eqn{\tilde\theta = \hat\theta_0 + \frac1n\sum_i \mathrm{IC}_i(\hat\theta_0)} :
#' un estimateur initial grossier corrigé de sa fonction d'influence devient
#' **asymptotiquement efficace**. C'est la forme abstraite de la correction du
#' lasso débiaisé (Module 22) et du score orthogonal de Neyman (Module 16).
#'
#' @param theta0 estimateur initial (scalaire ou vecteur).
#' @param ic_at fonction `theta -> IC(theta)` : vecteur de longueur n (theta
#'   scalaire) ou matrice n x p (theta vectoriel).
#' @return l'estimateur corrigé \eqn{\tilde\theta}.
#' @export
onestep <- function(theta0, ic_at) {
  ic <- ic_at(theta0)
  theta0 + (if (is.matrix(ic)) colMeans(ic) else mean(ic))
}

#' Fonction d'influence EFFICACE de l'ATE (doublement robuste / AIPW)
#'
#' \eqn{\psi_i = \mu_1(x_i)-\mu_0(x_i) + \frac{d_i(y_i-\mu_1)}{e_i}
#' - \frac{(1-d_i)(y_i-\mu_0)}{1-e_i} - \tau}. Son second moment est la **borne
#' d'efficacité semiparamétrique** de l'ATE ; l'AIPW/DML (Module 16) l'atteint.
#'
#' @param y résultat ; @param d traitement (0/1).
#' @param mu1,mu0 espérances conditionnelles estimées \eqn{E[Y|X,D=1/0]}.
#' @param e score de propension estimé \eqn{P(D=1|X)}.
#' @return liste : `ate`, `eif`, `se` (= sqrt(borne d'efficacité / n)).
#' @export
eif_ate <- function(y, d, mu1, mu0, e) {
  e <- pmin(pmax(e, 1e-6), 1 - 1e-6)
  score <- mu1 - mu0 + d * (y - mu1) / e - (1 - d) * (y - mu0) / (1 - e)
  ate <- mean(score); eif <- score - ate
  list(ate = ate, eif = eif, se = sqrt(mean(eif^2) / length(y)))
}
