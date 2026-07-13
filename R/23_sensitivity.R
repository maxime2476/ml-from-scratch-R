# =============================================================================
# Module 23 — Analyse de sensibilité à la confusion inobservée (Cinelli-Hazlett)
# Implémente les équations de derivations/23_sensitivity.qmd. R base.
# =============================================================================

#' R² partiel d'un effet à partir de sa statistique t
#'
#' \eqn{R^2_{Y\sim D\mid X} = t^2/(t^2+\mathrm{df})}.
#'
#' @param t statistique t de l'effet.
#' @param df degrés de liberté résiduels.
#' @return le R² partiel.
#' @export
partial_r2 <- function(t, df) t^2 / (t^2 + df)

#' Robustness value (Cinelli-Hazlett 2020)
#'
#' Force minimale (en R² partiel) qu'un confondeur inobservé — supposé aussi
#' associé au traitement qu'au résultat — devrait avoir pour **réduire l'effet
#' de q×100 %** (`alpha = 1`), ou pour le rendre **non significatif au seuil
#' alpha** (`alpha < 1`).
#'
#' @param t statistique t de l'effet.
#' @param df degrés de liberté résiduels.
#' @param q fraction de réduction visée (1 = ramener à 0).
#' @param alpha seuil de significativité (1 = point seulement).
#' @return la robustness value dans [0,1].
#' @export
robustness_value <- function(t, df, q = 1, alpha = 1) {
  fq <- q * abs(t) / sqrt(df)
  f_crit <- abs(qt(alpha / 2, df = df - 1)) / sqrt(df - 1)
  fqa <- fq - f_crit
  if (fqa < 0) return(0)
  rv <- 0.5 * (sqrt(fqa^4 + 4 * fqa^2) - fqa^2)
  # borne haute : quand la force dépasse ce que le confondeur peut expliquer
  if (fqa > 0 && fq > 1 / f_crit && f_crit > 0) {
    rvx <- (fq^2 - f_crit^2) / (1 + fq^2); rv <- max(rv, rvx)
  }
  min(rv, 1)
}

#' Estimation ajustée pour un confondeur de force donnée (OVB)
#'
#' Décale l'effet et l'erreur standard selon les R² partiels **hypothétiques** du
#' confondeur avec le traitement (`r2dz`) et avec le résultat (`r2yz`).
#'
#' @param estimate estimation de l'effet.
#' @param se erreur standard.
#' @param df degrés de liberté résiduels.
#' @param r2dz R² partiel confondeur-traitement.
#' @param r2yz R² partiel confondeur-résultat.
#' @param reduce réduire l'effet vers 0 (défaut TRUE).
#' @return liste : `estimate`, `se`, `t`, `bias`.
#' @export
adjusted_estimate <- function(estimate, se, df, r2dz, r2yz, reduce = TRUE) {
  bias <- sqrt(r2yz * r2dz / (1 - r2dz)) * se * sqrt(df)
  adj <- if (reduce) sign(estimate) * (abs(estimate) - bias) else sign(estimate) * (abs(estimate) + bias)
  adj_se <- se * sqrt((1 - r2yz) / (1 - r2dz)) * sqrt(df / (df - 1))
  list(estimate = adj, se = adj_se, t = adj / adj_se, bias = bias)
}

#' Analyse de sensibilité complète d'un effet OLS
#'
#' @param fit objet `ols` (Module 1).
#' @param treatment nom du régresseur d'intérêt.
#' @param q fraction de réduction pour la robustness value.
#' @param alpha seuil pour la robustness value significative.
#' @return liste : `estimate`, `se`, `t`, `df`, `r2yd`, `rv_q`, `rv_qa`.
#' @export
sensitivity_ols <- function(fit, treatment, q = 1, alpha = 0.05) {
  sm <- ols_summary(fit)$coefficients
  est <- sm[treatment, "estimate"]; se <- sm[treatment, "se"]; t <- sm[treatment, "t"]
  df <- fit$df.residual
  list(estimate = est, se = se, t = t, df = df,
       r2yd = partial_r2(t, df),
       rv_q = robustness_value(t, df, q = q, alpha = 1),
       rv_qa = robustness_value(t, df, q = q, alpha = alpha))
}
