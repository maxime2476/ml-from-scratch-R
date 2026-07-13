# =============================================================================
# Module 4 — Régularisation : ridge et lasso
# Implémente les équations de derivations/04_regularisation.qmd.
# Le coordinate descent + soft-thresholding du lasso est RÉIMPLÉMENTÉ localement
# (parti pris d'auto-suffisance : pas d'appel à optim_cd du Module 0).
# =============================================================================

# ---- Standardisation (écart-type de population, comme glmnet/lm.ridge) -------
.standardize <- function(X, center = TRUE, scale = TRUE) {
  X <- as.matrix(X); n <- nrow(X)
  ctr <- if (center) colMeans(X) else rep(0, ncol(X))
  Xc <- sweep(X, 2, ctr, "-")
  scl <- if (scale) sqrt(colSums(Xc^2) / n) else rep(1, ncol(X))  # sd de population
  scl[scl == 0] <- 1
  list(Xs = sweep(Xc, 2, scl, "/"), center = ctr, scale = scl)
}

#' Régression ridge (forme fermée, éq. 4.2)
#'
#' Résout \eqn{(X^TX+\lambda I)\beta = X^Ty} sur données éventuellement
#' standardisées, puis re-transforme les coefficients sur l'échelle d'origine.
#' L'intercept n'est pas pénalisé.
#'
#' @param X matrice de design n x p (sans colonne de constante).
#' @param y réponse.
#' @param lambda pénalité \eqn{\lambda \ge 0}.
#' @param standardize centrer-réduire les colonnes de X (défaut TRUE).
#' @param intercept inclure un intercept non pénalisé (défaut TRUE).
#' @return liste : `coefficients` (avec intercept si demandé), `beta` (pentes),
#'   `intercept`, `lambda`, `fitted`.
#' @export
ridge_fit <- function(X, y, lambda, standardize = TRUE, intercept = TRUE) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  st <- .standardize(X, center = intercept, scale = standardize)
  Xs <- st$Xs
  ybar <- if (intercept) mean(y) else 0
  yc <- y - ybar
  A <- crossprod(Xs) + lambda * diag(p)          # X'X + lambda I
  beta_s <- as.numeric(solve(A, crossprod(Xs, yc)))
  beta <- beta_s / st$scale                        # retour à l'échelle d'origine
  names(beta) <- colnames(X)
  b0 <- if (intercept) ybar - sum(st$center * beta) else 0
  fitted <- as.numeric(X %*% beta) + b0
  coefs <- if (intercept) c(`(Intercept)` = b0, beta) else beta
  list(coefficients = coefs, beta = beta, intercept = b0,
       lambda = lambda, fitted = fitted)
}

#' Biais, variance et EQM analytiques du ridge via la SVD (éq. 4.4-4.5)
#'
#' Décompose l'EQM du ridge en biais² + variance dans la base des composantes
#' principales (rotation V de la SVD de X). Sert à confronter la théorie au
#' Monte Carlo.
#'
#' @param X design (standardisé si l'on veut la cohérence avec `ridge_fit`).
#' @param beta_true vecteur des coefficients vrais (échelle des colonnes de X).
#' @param sigma2 variance des erreurs.
#' @param lambda pénalité.
#' @return liste : `bias2`, `variance`, `mse` (totaux), et vecteurs par composante.
#' @export
ridge_bias_var <- function(X, beta_true, sigma2, lambda) {
  sv <- svd(X)
  d <- sv$d
  alpha <- as.numeric(crossprod(sv$v, beta_true))   # alpha = V' beta
  f <- d^2 / (d^2 + lambda)                           # rétrécissement (éq. 4.3)
  bias_alpha <- (f - 1) * alpha                       # éq. 4.4
  var_alpha  <- sigma2 * d^2 / (d^2 + lambda)^2        # éq. 4.4
  list(bias2 = sum(bias_alpha^2), variance = sum(var_alpha),
       mse = sum(bias_alpha^2 + var_alpha),
       bias_alpha = bias_alpha, var_alpha = var_alpha, d = d, shrink = f)
}

#' Opérateur de soft-thresholding (éq. 4.9)
#'
#' \eqn{\mathcal S(z,\lambda) = \mathrm{sign}(z)\,(|z|-\lambda)_+}. Vectorisé.
#'
#' @param z scalaire ou vecteur.
#' @param lambda seuil (>= 0).
#' @return la (les) valeur(s) seuillée(s).
#' @export
soft_threshold <- function(z, lambda) sign(z) * pmax(abs(z) - lambda, 0)

#' Lasso par coordinate descent (éq. 4.10)
#'
#' Boucle de descente par coordonnées avec mise à jour du résidu ; chaque
#' coordonnée est mise à jour par soft-thresholding (Prop. 4.2). L'intercept
#' n'est pas pénalisé.
#'
#' @param X matrice de design n x p (sans constante).
#' @param y réponse.
#' @param lambda pénalité \eqn{\lambda \ge 0} (convention \eqn{\tfrac12\|y-X\beta\|^2 + \lambda\|\beta\|_1}).
#' @param standardize centrer-réduire X (défaut TRUE).
#' @param intercept intercept non pénalisé (défaut TRUE).
#' @param maxit balayages maximaux.
#' @param tol tolérance d'arrêt (variation max des coefficients).
#' @return liste : `coefficients`, `beta`, `intercept`, `lambda`, `iter`, `fitted`.
#' @export
lasso_fit <- function(X, y, lambda, standardize = TRUE, intercept = TRUE,
                      maxit = 1e4L, tol = 1e-9) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  st <- .standardize(X, center = intercept, scale = standardize)
  Xs <- st$Xs
  ybar <- if (intercept) mean(y) else 0
  yc <- y - ybar

  v <- colSums(Xs^2)                               # ||x_j||^2
  beta <- rep(0, p)
  r <- yc                                           # résidu = yc - Xs beta (beta=0)
  it <- 0L
  for (it in seq_len(maxit)) {
    max_change <- 0
    for (j in seq_len(p)) {
      zj <- sum(Xs[, j] * r) + v[j] * beta[j]       # x_j' r^{(j)} (éq. 4.10)
      bj <- soft_threshold(zj, lambda) / v[j]        # Prop. 4.2
      d  <- bj - beta[j]
      if (d != 0) { r <- r - Xs[, j] * d; beta[j] <- bj
                    max_change <- max(max_change, abs(d)) }
    }
    if (max_change < tol) break
  }
  beta_orig <- beta / st$scale
  names(beta_orig) <- colnames(X)
  b0 <- if (intercept) ybar - sum(st$center * beta_orig) else 0
  fitted <- as.numeric(X %*% beta_orig) + b0
  coefs <- if (intercept) c(`(Intercept)` = b0, beta_orig) else beta_orig
  list(coefficients = coefs, beta = beta_orig, intercept = b0,
       lambda = lambda, iter = it, fitted = fitted)
}
