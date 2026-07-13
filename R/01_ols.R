# =============================================================================
# Module 1 — Moindres carrés ordinaires et inférence
# Implémente les équations de derivations/01_ols.qmd.
#
# S'appuie sur la QR de Householder du Module 0 (solve_ls_qr) : les
# factorisations d'algèbre linéaire sont les primitives numériques partagées du
# projet. Interface "formule" façon lm().
# =============================================================================

#' Ajustement MCO par QR, avec inférence
#'
#' Résout les équations normales (éq. 1.2) via la QR du Module 0, puis calcule la
#' matrice de variance classique \eqn{s^2(X^TX)^{-1}} (éq. 1.5), l'estimateur
#' sans biais \eqn{s^2} (éq. 1.6) et les leviers (diagonale de la matrice
#' chapeau). Aucune inversion de \eqn{X^TX} n'est formée : on exploite le facteur
#' \eqn{R} de la QR (\eqn{X^TX = R^T R}, donc \eqn{(X^TX)^{-1}=R^{-1}R^{-T}}).
#'
#' @param formula formule façon `lm` (ex. `y ~ x1 + x2`).
#' @param data data.frame contenant les variables de la formule.
#' @return objet de classe `ols` : `coefficients`, `vcov`, `sigma2`,
#'   `df.residual`, `fitted`, `residuals`, `hat`, `rss`, `tss`, `rank`, `n`,
#'   `p`, `has_intercept`, `terms`, `xlevels`.
#' @export
ols_fit <- function(formula, data) {
  mf <- model.frame(formula, data)
  tt <- attr(mf, "terms")
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(tt, mf)
  n <- nrow(X); p <- ncol(X)
  if (n <= p) stop("n <= p : degrés de liberté résiduels insuffisants.")

  fit <- solve_ls_qr(X, y)                 # éq. (1.3) via QR (Module 0)
  R1  <- fit$R                              # bloc triangulaire p x p, X'X = R1'R1
  beta <- fit$coefficients
  names(beta) <- colnames(X)

  # (X'X)^{-1} = R1^{-1} R1^{-T}, sans former X'X ---------------------------
  Rinv  <- backsolve(R1, diag(p))          # inverse d'une triangulaire supérieure
  XtXinv <- Rinv %*% t(Rinv)
  dimnames(XtXinv) <- list(colnames(X), colnames(X))

  df.residual <- n - p
  sigma2 <- fit$rss / df.residual          # éq. (1.6)
  vcov   <- sigma2 * XtXinv                 # éq. (1.5)

  # Leviers : diag(H) = rowSums((X Rinv)^2) car H = (X R1^{-1})(X R1^{-1})^T ---
  W   <- X %*% Rinv
  hat <- rowSums(W^2)

  has_intercept <- attr(tt, "intercept") == 1
  ybar <- if (has_intercept) mean(y) else 0
  tss  <- sum((y - ybar)^2)

  structure(list(
    coefficients = beta, vcov = vcov, sigma2 = sigma2, sigma = sqrt(sigma2),
    df.residual = df.residual, fitted = fit$fitted, residuals = fit$residuals,
    hat = hat, rss = fit$rss, tss = tss, rank = p, n = n, p = p,
    has_intercept = has_intercept, terms = tt, XtXinv = XtXinv,
    model_matrix = X, response = y
  ), class = "ols")
}

#' Tableau récapitulatif d'un ajustement MCO
#'
#' Statistiques t et p-values (éq. 1.8), \eqn{R^2} et \eqn{\bar R^2}
#' (éq. 1.11-1.12), et test F global (éq. 1.9, cas pente nulle). Reproduit
#' `summary.lm()`.
#'
#' @param object objet `ols`.
#' @return liste : `coefficients` (data.frame estimate/se/t/p_value), `r2`,
#'   `adj_r2`, `sigma`, `fstatistic` (value, numdf, dendf, p_value).
#' @export
ols_summary <- function(object) {
  b  <- object$coefficients
  se <- sqrt(diag(object$vcov))
  tval <- b / se
  pval <- 2 * pt(-abs(tval), df = object$df.residual)   # éq. (1.8)
  coefs <- data.frame(estimate = b, se = se, t = tval, p_value = pval,
                      row.names = names(b))

  r2 <- if (object$has_intercept) 1 - object$rss / object$tss else NA_real_
  adj_r2 <- if (object$has_intercept)
    1 - (1 - r2) * (object$n - 1) / object$df.residual else NA_real_   # éq. (1.12)

  # F global : toutes les pentes nulles (df1 = p - has_intercept) -----------
  fstat <- NULL
  if (object$has_intercept && object$p > 1) {
    df1 <- object$p - 1
    Fval <- (r2 / df1) / ((1 - r2) / object$df.residual)
    fstat <- list(value = Fval, numdf = df1, dendf = object$df.residual,
                  p_value = pf(Fval, df1, object$df.residual, lower.tail = FALSE))
  }
  list(coefficients = coefs, r2 = r2, adj_r2 = adj_r2, sigma = object$sigma,
       df.residual = object$df.residual, fstatistic = fstat)
}

#' Intervalles de confiance des coefficients (éq. 1.8)
#'
#' @param object objet `ols`.
#' @param level niveau de confiance (défaut 0.95).
#' @return matrice à deux colonnes (bornes inf/sup).
#' @export
ols_confint <- function(object, level = 0.95) {
  b  <- object$coefficients
  se <- sqrt(diag(object$vcov))
  q  <- qt(1 - (1 - level) / 2, df = object$df.residual)
  ci <- cbind(b - q * se, b + q * se)
  colnames(ci) <- paste0(c((1 - level) / 2, 1 - (1 - level) / 2) * 100, " %")
  rownames(ci) <- names(b)
  ci
}

#' Test F de restrictions linéaires R beta = r (éq. 1.9)
#'
#' @param object objet `ols`.
#' @param R matrice q x p des restrictions (rang q).
#' @param r vecteur cible (longueur q ; défaut zéro).
#' @return liste : `F`, `q`, `df.residual`, `p_value`.
#' @export
ols_ftest <- function(object, R, r = rep(0, nrow(R))) {
  R <- as.matrix(R)
  if (ncol(R) != object$p) stop("R doit avoir p colonnes.")
  q <- nrow(R)
  Rb_r <- as.numeric(R %*% object$coefficients - r)
  mid  <- R %*% object$XtXinv %*% t(R)              # R (X'X)^{-1} R'
  Fval <- as.numeric(t(Rb_r) %*% solve(mid, Rb_r)) / (q * object$sigma2)
  list(F = Fval, q = q, df.residual = object$df.residual,
       p_value = pf(Fval, q, object$df.residual, lower.tail = FALSE))
}

#' Coefficient FWL : régression de M1 y sur M1 X2 (éq. 1.10)
#'
#' Renvoie \eqn{\hat\beta_2} obtenu en résidualisant y et X2 par rapport à X1
#' (théorème de Frisch-Waugh-Lovell). Sert à la vérification numérique : ce
#' vecteur doit égaler le sous-vecteur \eqn{\hat\beta_2} de la régression complète.
#'
#' @param y réponse (vecteur).
#' @param X1 bloc à partialiser (matrice n x p1).
#' @param X2 bloc d'intérêt (matrice n x p2).
#' @return le vecteur \eqn{\hat\beta_2}.
#' @export
fwl_beta2 <- function(y, X1, X2) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2); y <- as.numeric(y)
  # Résidualisation par MCO sur X1 (via QR du Module 0)
  resid_on_X1 <- function(z) z - X1 %*% solve_ls_qr(X1, z)$coefficients
  M1y  <- resid_on_X1(y)
  M1X2 <- apply(X2, 2, resid_on_X1)
  as.numeric(solve_ls_qr(M1X2, M1y)$coefficients)
}
