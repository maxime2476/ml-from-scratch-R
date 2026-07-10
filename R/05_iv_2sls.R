# =============================================================================
# Module 5 — Variables instrumentales et doubles moindres carrés (2SLS)
# Implémente les équations de derivations/05_iv_2sls.qmd.
# Interface matricielle : X = régresseurs (incl. constante, exogènes, endogènes),
# Z = instruments (incl. constante, exogènes, instruments exclus).
# Les projections passent par la QR du Module 0 (solve_ls_qr).
# =============================================================================

#' Estimateur 2SLS (forme fermée, éq. 5.4)
#'
#' Projette X sur l'espace des instruments (\eqn{\hat X = P_Z X}, étape 1) puis
#' régresse y sur \eqn{\hat X} (étape 2), d'où
#' \eqn{\hat\beta = (X^T P_Z X)^{-1} X^T P_Z y}. La variance (éq. 5.5) utilise les
#' résidus sur le VRAI X (équation structurelle), pas \eqn{\hat X}.
#'
#' @param y réponse (vecteur longueur n).
#' @param X matrice n x k des régresseurs (constante incluse ; colonnes
#'   exogènes et endogènes).
#' @param Z matrice n x m des instruments (m >= k ; inclut la constante et les
#'   régresseurs exogènes, plus les instruments exclus).
#' @return liste : `coefficients`, `vcov`, `se`, `sigma2`, `df.residual`,
#'   `fitted`, `residuals`, `Xhat`.
tsls_fit <- function(y, X, Z) {
  X <- as.matrix(X); Z <- as.matrix(Z); y <- as.numeric(y)
  n <- nrow(X); k <- ncol(X)
  if (nrow(Z) != n) stop("Z et X doivent avoir n lignes.")
  if (ncol(Z) < k) stop("Condition d'ordre violée : ncol(Z) < ncol(X).")

  # Étape 1 : Xhat = P_Z X (valeurs ajustées de X régressé sur Z, colonne par colonne)
  Xhat <- apply(X, 2, function(col) solve_ls_qr(Z, col)$fitted)
  Xhat <- matrix(Xhat, n, k, dimnames = list(NULL, colnames(X)))

  # Étape 2 : régression de y sur Xhat -> beta = (Xhat'Xhat)^{-1} Xhat'y
  fit2 <- solve_ls_qr(Xhat, y)
  beta <- fit2$coefficients; names(beta) <- colnames(X)

  # Variance : résidus sur le VRAI X (éq. 5.5)
  fitted <- as.numeric(X %*% beta)
  resid  <- y - fitted
  df.residual <- n - k
  sigma2 <- sum(resid^2) / df.residual
  Rinv <- backsolve(fit2$R, diag(k))
  XtPzXinv <- Rinv %*% t(Rinv)                 # (Xhat'Xhat)^{-1} = (X'P_Z X)^{-1}
  vcov <- sigma2 * XtPzXinv
  dimnames(vcov) <- list(colnames(X), colnames(X))

  list(coefficients = beta, vcov = vcov, se = sqrt(diag(vcov)), sigma2 = sigma2,
       df.residual = df.residual, fitted = fitted, residuals = resid, Xhat = Xhat)
}

#' Statistique F de première étape (force des instruments)
#'
#' Teste la nullité conjointe des coefficients des instruments EXCLUS dans la
#' régression de première étape du régresseur endogène sur Z (diagnostic
#' d'instruments faibles ; règle empirique F < 10, cf. §5 de la dérivation).
#'
#' @param x_endog régresseur endogène (vecteur).
#' @param Z matrice des instruments (constante + exogènes + instruments exclus).
#' @param excluded indices des colonnes de Z correspondant aux instruments exclus.
#' @return liste : `F`, `df1`, `df2`, `p_value`.
first_stage_F <- function(x_endog, Z, excluded) {
  Z <- as.matrix(Z); x <- as.numeric(x_endog)
  n <- nrow(Z); m <- ncol(Z); q <- length(excluded)
  rss_full <- solve_ls_qr(Z, x)$rss                       # première étape complète
  rss_red  <- solve_ls_qr(Z[, -excluded, drop = FALSE], x)$rss  # sans instruments exclus
  Fval <- ((rss_red - rss_full) / q) / (rss_full / (n - m))
  list(F = Fval, df1 = q, df2 = n - m,
       p_value = pf(Fval, q, n - m, lower.tail = FALSE))
}
