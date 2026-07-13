# =============================================================================
# Module 21 — Données de panel et effets fixes
# Implémente les équations de derivations/21_panel.qmd. R base + Module 0.
# =============================================================================

#' Estimateur à effets fixes (within) pour données de panel (éq. 21.2)
#'
#' Centre les variables par unité (transformation within, qui élimine l'effet
#' fixe \eqn{\alpha_i}) puis applique l'OLS. Renvoie aussi la variance groupée
#' (clustered) par unité (éq. 21.3). Équivaut au LSDV (Prop. 21.1).
#'
#' @param formula formule des régresseurs variables (sans les indicatrices).
#' @param data data.frame.
#' @param id nom de la colonne identifiant l'unité.
#' @return liste : `coefficients`, `vcov` (classique within), `vcov_cluster`,
#'   `se`, `se_cluster`, `df.residual`, `residuals`, `N`, `NT`.
#' @export
fe_within <- function(formula, data, id) {
  mf <- model.frame(formula, data)
  tt <- attr(mf, "terms")
  y <- as.numeric(model.response(mf))
  X <- model.matrix(tt, mf)
  # retirer l'intercept (absorbé par les effets fixes)
  if ("(Intercept)" %in% colnames(X)) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  g <- as.factor(data[[id]][as.integer(rownames(mf))])
  N <- nlevels(g); NT <- nrow(X); k <- ncol(X)

  demean <- function(M) M - rowsum(M, g)[as.integer(g), , drop = FALSE] / as.numeric(table(g))[as.integer(g)]
  Xt <- demean(X); yt <- demean(matrix(y, ncol = 1))[, 1]

  XtX <- crossprod(Xt); XtXinv <- solve(XtX)
  beta <- as.numeric(XtXinv %*% crossprod(Xt, yt)); names(beta) <- colnames(X)
  resid <- yt - as.numeric(Xt %*% beta)
  df <- NT - N - k
  sigma2 <- sum(resid^2) / df
  vcov <- sigma2 * XtXinv                                # within classique
  dimnames(vcov) <- list(colnames(X), colnames(X))

  # Viande groupée par unité (éq. 21.3), avec correction d.d.l. façon sandwich/plm.
  meat <- matrix(0, k, k)
  for (lev in levels(g)) {
    ix <- which(g == lev)
    Xi <- Xt[ix, , drop = FALSE]; ui <- resid[ix]
    s <- crossprod(Xi, ui)
    meat <- meat + s %*% t(s)
  }
  adj <- (N / (N - 1)) * ((NT - 1) / (NT - k))          # correction de vcovCL (HC1-like)
  vcov_cl <- adj * (XtXinv %*% meat %*% XtXinv)
  dimnames(vcov_cl) <- list(colnames(X), colnames(X))

  list(coefficients = beta, vcov = vcov, vcov_cluster = vcov_cl,
       se = sqrt(diag(vcov)), se_cluster = sqrt(diag(vcov_cl)),
       df.residual = df, residuals = resid, N = N, NT = NT)
}
