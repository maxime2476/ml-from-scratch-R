# =============================================================================
# Module 2 — Hétéroscédasticité, erreurs robustes, GLS/WLS
# Implémente les équations de derivations/02_gls_robust.qmd.
# Les fonctions HC/NW prennent un objet `ols` (Module 1) ; wls_fit/gls_fit ont
# une interface formule et réutilisent la QR du Module 0.
# =============================================================================

#' Matrice de variance sandwich générique (éq. 2.1)
#'
#' Calcule \eqn{(X^TX)^{-1}(X^T \Omega X)(X^TX)^{-1}} pour une « viande »
#' \eqn{X^T \Omega X} fournie. Sert de brique commune à HC et Newey-West.
#'
#' @param fit objet `ols`.
#' @param meat matrice p x p (la viande \eqn{X^T \Omega X} estimée).
#' @return matrice de variance p x p.
sandwich_vcov <- function(fit, meat) {
  bread <- fit$XtXinv
  V <- bread %*% meat %*% bread
  dimnames(V) <- dimnames(fit$vcov)
  V
}

#' Variances robustes à l'hétéroscédasticité HC0–HC3 (éq. 2.2)
#'
#' Viande \eqn{\sum_i \phi_i \hat\varepsilon_i^2 x_i x_i^T} avec le facteur
#' \eqn{\phi_i} dépendant du type (Prop. 2.2 pour la correction de levier).
#'
#' @param fit objet `ols`.
#' @param type "HC0", "HC1", "HC2" ou "HC3".
#' @return matrice de variance p x p (compatible `sandwich::vcovHC`).
vcov_hc <- function(fit, type = c("HC3", "HC0", "HC1", "HC2")) {
  type <- match.arg(type)
  X <- fit$model_matrix
  e <- fit$residuals
  h <- fit$hat
  n <- fit$n; p <- fit$p
  phi <- switch(type,
    HC0 = rep(1, n),
    HC1 = rep(n / (n - p), n),
    HC2 = 1 / (1 - h),
    HC3 = 1 / (1 - h)^2)
  omega <- phi * e^2                       # \phi_i \hat\varepsilon_i^2
  meat  <- crossprod(X, omega * X)          # X^T diag(omega) X
  sandwich_vcov(fit, meat)
}

#' Tests t robustes (coeftest) avec une variance donnée
#'
#' Reproduit `lmtest::coeftest` : t = coef / se_robuste, loi de Student à
#' `df.residual` degrés de liberté.
#'
#' @param fit objet `ols`.
#' @param vcov matrice de variance (p.ex. sortie de `vcov_hc`).
#' @return data.frame estimate/se/t/p_value.
coeftest_hc <- function(fit, vcov = vcov_hc(fit, "HC3")) {
  b  <- fit$coefficients
  se <- sqrt(diag(vcov))
  tval <- b / se
  pval <- 2 * pt(-abs(tval), df = fit$df.residual)
  data.frame(estimate = b, se = se, t = tval, p_value = pval,
             row.names = names(b))
}

#' Variance HAC de Newey-West (éq. 2.6)
#'
#' Viande \eqn{\hat\Gamma_0 + \sum_{\ell=1}^L w_\ell(\hat\Gamma_\ell +
#' \hat\Gamma_\ell^T)} avec poids de Bartlett \eqn{w_\ell = 1 - \ell/(L+1)}.
#' Calé sur `sandwich::NeweyWest(..., prewhite = FALSE, adjust = FALSE)`.
#'
#' @param fit objet `ols`.
#' @param lag nombre de retards L.
#' @return matrice de variance p x p.
vcov_nw <- function(fit, lag) {
  X <- fit$model_matrix
  e <- fit$residuals
  n <- fit$n
  U <- X * e                                # lignes u_i^T = e_i x_i^T
  S <- crossprod(U)                          # \hat\Gamma_0 = sum u_i u_i^T
  for (l in seq_len(lag)) {
    w <- 1 - l / (lag + 1)
    G <- crossprod(U[(l + 1):n, , drop = FALSE], U[1:(n - l), , drop = FALSE])
    S <- S + w * (G + t(G))                  # \Gamma_l + \Gamma_l^T
  }
  sandwich_vcov(fit, S)
}

#' Moindres carrés pondérés (WLS, éq. 2.5)
#'
#' Résout \eqn{(X^T W X)^{-1} X^T W y} en appliquant l'OLS (QR, Module 0) aux
#' données transformées \eqn{\sqrt{w_i}\,(x_i, y_i)}. Reproduit `lm(weights=)`.
#'
#' @param formula formule façon `lm`.
#' @param data data.frame.
#' @param weights vecteur de poids \eqn{w_i = 1/\sigma_i^2} (longueur n).
#' @return objet de classe `wls` (mêmes champs qu'`ols`).
wls_fit <- function(formula, data, weights) {
  mf <- model.frame(formula, data)
  tt <- attr(mf, "terms")
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(tt, mf)
  n <- nrow(X); p <- ncol(X)
  w <- as.numeric(weights)
  if (length(w) != n) stop("length(weights) doit valoir n.")
  sw <- sqrt(w)
  fit <- solve_ls_qr(X * sw, y * sw)        # OLS sur données transformées
  beta <- fit$coefficients; names(beta) <- colnames(X)

  resid <- y - as.numeric(X %*% beta)        # résidus non pondérés (comme lm)
  wrss  <- sum(w * resid^2)
  df.residual <- n - p
  sigma2 <- wrss / df.residual               # sigma^2 pondéré
  Rinv <- backsolve(fit$R, diag(p))
  XtWXinv <- Rinv %*% t(Rinv)                 # (X^T W X)^{-1}
  dimnames(XtWXinv) <- list(colnames(X), colnames(X))
  vcov <- sigma2 * XtWXinv

  structure(list(
    coefficients = beta, vcov = vcov, sigma2 = sigma2, sigma = sqrt(sigma2),
    df.residual = df.residual, fitted = as.numeric(X %*% beta), residuals = resid,
    weights = w, rss = wrss, n = n, p = p, XtXinv = XtWXinv,
    model_matrix = X, response = y
  ), class = c("wls", "ols"))
}

#' Moindres carrés généralisés (GLS, éq. 2.4) avec Omega connue
#'
#' Résout \eqn{(X^T \Omega^{-1} X)^{-1} X^T \Omega^{-1} y} en appliquant l'OLS aux
#' données transformées \eqn{P^{-1}(X, y)} où \eqn{\Omega = P P^T} (Cholesky,
#' Module 0). Théorème d'Aitken (Th. 2.3).
#'
#' @param formula formule façon `lm`.
#' @param data data.frame.
#' @param Omega matrice n x n SPD (structure de covariance des erreurs, à un
#'   facteur d'échelle sigma^2 près).
#' @return objet de classe `gls_fit` (mêmes champs qu'`ols`).
gls_fit <- function(formula, data, Omega) {
  mf <- model.frame(formula, data)
  tt <- attr(mf, "terms")
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(tt, mf)
  n <- nrow(X); p <- ncol(X)
  if (!all(dim(Omega) == n)) stop("Omega doit être n x n.")

  P  <- chol_crout(Omega)                    # Omega = P P^T (P triangulaire inf.)
  # Transformation P^{-1} par descente triangulaire (colonne par colonne).
  Pinv_mult <- function(Z) apply(as.matrix(Z), 2, function(col) forward_substitution(P, col))
  Xs <- Pinv_mult(X)                          # P^{-1} X
  ys <- as.numeric(forward_substitution(P, y))# P^{-1} y

  fit <- solve_ls_qr(Xs, ys)
  beta <- fit$coefficients; names(beta) <- colnames(X)
  df.residual <- n - p
  # sigma^2 = (transformed RSS)/(n-p) = resid' Omega^{-1} resid /(n-p)
  sigma2 <- fit$rss / df.residual
  Rinv <- backsolve(fit$R, diag(p))
  XtOinvXinv <- Rinv %*% t(Rinv)              # (X^T Omega^{-1} X)^{-1}
  dimnames(XtOinvXinv) <- list(colnames(X), colnames(X))
  vcov <- sigma2 * XtOinvXinv

  structure(list(
    coefficients = beta, vcov = vcov, sigma2 = sigma2, sigma = sqrt(sigma2),
    df.residual = df.residual,
    fitted = as.numeric(X %*% beta), residuals = y - as.numeric(X %*% beta),
    n = n, p = p, XtXinv = XtOinvXinv, model_matrix = X, response = y
  ), class = c("gls_fit", "ols"))
}
