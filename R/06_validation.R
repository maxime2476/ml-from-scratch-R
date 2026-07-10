# =============================================================================
# Module 6 — Validation de modèles
# Implémente les équations de derivations/06_validation.qmd.
# S'appuie sur la QR du Module 0 pour les ajustements linéaires.
# =============================================================================

#' LOOCV fermé pour la régression linéaire (éq. 6.3)
#'
#' Calcule \eqn{\mathrm{CV}_n = \frac1n\sum_i (\hat\varepsilon_i/(1-h_{ii}))^2} en
#' UN seul ajustement, via les leviers \eqn{h_{ii}} (Théorème 6.2).
#'
#' @param X matrice de design (constante incluse).
#' @param y réponse.
#' @return liste : `cv` (erreur LOOCV moyenne), `loo_resid`, `h` (leviers).
loocv_linear <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y); p <- ncol(X)
  fit <- solve_ls_qr(X, y)
  Rinv <- backsolve(fit$R, diag(p))
  h <- rowSums((X %*% Rinv)^2)              # diag(H)
  loo <- fit$residuals / (1 - h)
  list(cv = mean(loo^2), loo_resid = loo, h = h)
}

#' Validation croisée généralisée (GCV)
#'
#' Remplace \eqn{h_{ii}} par la moyenne \eqn{\bar h = p/n} dans la formule LOOCV.
#'
#' @param X matrice de design (constante incluse).
#' @param y réponse.
#' @return liste : `gcv`.
gcv_linear <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  fit <- solve_ls_qr(X, y)
  list(gcv = mean((fit$residuals / (1 - p / n))^2))
}

#' Validation croisée K-fold générique (éq. 6.2)
#'
#' @param X matrice de design (constante incluse si le modèle en a une).
#' @param y réponse.
#' @param K nombre de blocs (défaut 10).
#' @param fit_fun fonction `(Xtr, ytr) -> modèle` (défaut : OLS via QR).
#' @param pred_fun fonction `(modèle, Xte) -> prédictions` (défaut : OLS).
#' @param seed graine pour le tirage des blocs.
#' @return liste : `cv`, `se`, `fold_errors`, `folds`.
kfold_cv <- function(X, y, K = 10L,
                     fit_fun = function(Xtr, ytr) solve_ls_qr(Xtr, ytr)$coefficients,
                     pred_fun = function(beta, Xte) as.numeric(Xte %*% beta),
                     seed = NULL) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X)
  if (!is.null(seed)) set.seed(seed)
  folds <- sample(rep_len(seq_len(K), n))
  err <- numeric(n)
  for (k in seq_len(K)) {
    te <- which(folds == k); tr <- which(folds != k)
    m <- fit_fun(X[tr, , drop = FALSE], y[tr])
    err[te] <- (y[te] - pred_fun(m, X[te, , drop = FALSE]))^2
  }
  list(cv = mean(err), se = sd(err) / sqrt(n), fold_errors = err, folds = folds)
}

#' Cp de Mallows (éq. 6.4)
#'
#' \eqn{C_p = \overline{\mathrm{err}} + 2\sigma^2 p/n}, correction de l'optimisme.
#'
#' @param fit objet `ols` (Module 1).
#' @param sigma2 estimation de la variance du bruit (p.ex. d'un modèle riche).
#' @return la valeur du Cp (échelle erreur quadratique moyenne).
mallows_cp <- function(fit, sigma2) {
  mean(fit$residuals^2) + 2 * sigma2 * fit$p / fit$n
}

# ---- Critères d'information (AIC/BIC) ----------------------------------------

#' Log-vraisemblance gaussienne (régression linéaire) au MLE
#' @param rss somme des carrés des résidus.
#' @param n taille d'échantillon.
#' @return la log-vraisemblance maximisée.
gaussian_loglik <- function(rss, n) -0.5 * n * (log(2 * pi) + log(rss / n) + 1)

#' AIC et BIC d'un modèle ajusté (éq. 6.5, 6.7)
#'
#' Gère les objets `ols` (Module 1 ; variance comptée comme paramètre, \eqn{k=p+1})
#' et `glm_irls` (Module 3 ; \eqn{k=p}). Reproduit `AIC()`/`BIC()` de R.
#'
#' @param fit objet `ols` ou `glm_irls`.
#' @return liste : `aic`, `bic`, `loglik`, `k`, `n`.
info_criteria <- function(fit) {
  if (inherits(fit, "glm_irls")) {
    ll <- fit[["loglik"]]; k <- fit[["rank"]]; n <- fit[["n"]]
  } else if (inherits(fit, "ols")) {
    n <- fit[["n"]]; k <- fit[["p"]] + 1L        # +1 pour sigma^2 (comme lm)
    ll <- gaussian_loglik(fit[["rss"]], n)
  } else stop("fit doit être de classe 'ols' ou 'glm_irls'.")
  list(aic = -2 * ll + 2 * k, bic = -2 * ll + log(n) * k,
       loglik = ll, k = k, n = n)
}

#' Estimation Monte Carlo de la décomposition biais-variance (éq. 6.1)
#'
#' Sur un DGP connu, estime irréductible / biais² / variance de la prédiction en
#' un point de test x0, pour un ajusteur donné.
#'
#' @param gen_y fonction `() -> y` générant une réponse (design X0 fixe, bruit aléatoire).
#' @param fit_pred fonction `(y) -> prédiction en x0` (scalaire).
#' @param f0 vraie valeur \eqn{f(x_0)}.
#' @param sigma2 variance du bruit (partie irréductible).
#' @param R nombre de réplications.
#' @return liste : `irreducible`, `bias2`, `variance`, `mse_pred`, `total`.
bias_variance_mc <- function(gen_y, fit_pred, f0, sigma2, R = 5000L) {
  preds <- numeric(R)
  for (r in seq_len(R)) preds[r] <- fit_pred(gen_y())
  fbar <- mean(preds)
  list(irreducible = sigma2, bias2 = (f0 - fbar)^2, variance = var(preds),
       total = sigma2 + (f0 - fbar)^2 + var(preds))
}
