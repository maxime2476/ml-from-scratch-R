# =============================================================================
# Module 19 — Prédiction conforme (split conformal)
# Implémente les équations de derivations/19_conformal.qmd. R base.
# =============================================================================

#' Quantile conforme (éq. 19.1)
#'
#' Renvoie le \eqn{\lceil(1-\alpha)(n+1)\rceil}-ème plus petit score
#' (\eqn{+\infty} si ce rang dépasse n).
#'
#' @param scores scores de non-conformité de calibration.
#' @param alpha niveau (couverture visée \eqn{1-\alpha}).
#' @return le quantile conforme \eqn{\hat q}.
#' @export
conformal_quantile <- function(scores, alpha = 0.1) {
  n <- length(scores)
  k <- ceiling((1 - alpha) * (n + 1))
  if (k > n) return(Inf)
  sort(scores)[k]
}

#' Prédiction conforme par découpage (split conformal)
#'
#' Ajuste un modèle sur l'échantillon d'entraînement, calibre les scores de
#' non-conformité \eqn{|y-\hat\mu(x)|} sur l'échantillon de calibration, et
#' renvoie des intervalles \eqn{\hat\mu(x)\pm\hat q} pour `X_test` (Théorème 19.1).
#' Avec `normalize = TRUE`, les scores sont divisés par une dispersion locale
#' \eqn{\hat\sigma(x)} (intervalles de largeur variable).
#'
#' @param X_train,y_train données d'entraînement.
#' @param X_cal,y_cal données de calibration.
#' @param X_test points où prédire.
#' @param fit_fn fonction `(X, y) -> modèle`.
#' @param predict_fn fonction `(modèle, X) -> prédictions`.
#' @param alpha niveau (défaut 0.1 -> 90 %).
#' @param normalize normalisation locale (défaut FALSE).
#' @param sigma_fn (si normalize) fonction `(X) -> dispersion locale positive`.
#' @return liste : `lower`, `upper`, `pred`, `qhat`.
#' @export
conformal_split <- function(X_train, y_train, X_cal, y_cal, X_test,
                            fit_fn, predict_fn, alpha = 0.1,
                            normalize = FALSE, sigma_fn = NULL) {
  model <- fit_fn(X_train, y_train)
  mu_cal <- predict_fn(model, X_cal)
  scores <- abs(as.numeric(y_cal) - mu_cal)
  if (normalize) {
    if (is.null(sigma_fn)) stop("normalize = TRUE requiert sigma_fn.")
    s_cal <- pmax(sigma_fn(X_cal), 1e-8); scores <- scores / s_cal
  }
  qhat <- conformal_quantile(scores, alpha)
  pred <- predict_fn(model, X_test)
  half <- if (normalize) qhat * pmax(sigma_fn(X_test), 1e-8) else qhat
  list(lower = pred - half, upper = pred + half, pred = pred, qhat = qhat)
}
