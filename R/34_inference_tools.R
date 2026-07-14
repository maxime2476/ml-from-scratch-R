# =============================================================================
# Module 34 — Outils d'inference : methode delta et tests multiples
# Implemente les equations de derivations/34_inference_tools.qmd. R base.
# Deux utilitaires transversaux : (a) l'ecart-type d'une FONCTION d'un estimateur
# (methode delta), (b) le controle des faux positifs quand on teste BEAUCOUP
# d'hypotheses (Bonferroni, Benjamini-Hochberg).
# =============================================================================

#' Methode delta : variance asymptotique d'une fonction d'un estimateur
#'
#' Si \eqn{\sqrt n(\hat\theta-\theta)\to\mathcal N(0,\Sigma)}, alors pour une
#' fonction reguliere \eqn{g}, \eqn{\sqrt n(g(\hat\theta)-g(\theta))\to\mathcal N
#' (0,\nabla g^\top\Sigma\nabla g)} (developpement de Taylor au premier ordre).
#' Le gradient est calcule par differences finies centrees.
#'
#' @param theta estimateur \eqn{\hat\theta} (vecteur).
#' @param vcov matrice de covariance estimee \eqn{\hat\Sigma}.
#' @param g fonction `theta -> scalaire`.
#' @param level niveau de confiance.
#' @return liste : `estimate`, `se`, `ci`.
#' @export
delta_method <- function(theta, vcov, g, level = 0.95) {
  eps <- 1e-6; k <- length(theta)
  grad <- vapply(seq_len(k), function(j) {
    tp <- theta; tp[j] <- tp[j] + eps; tm <- theta; tm[j] <- tm[j] - eps
    (g(tp) - g(tm)) / (2 * eps)
  }, numeric(1))
  est <- as.numeric(g(theta)); se <- sqrt(as.numeric(t(grad) %*% vcov %*% grad))
  z <- qnorm(1 - (1 - level) / 2)
  list(estimate = est, se = se, ci = est + c(-1, 1) * z * se)
}

#' Correction de Bonferroni (controle du FWER)
#'
#' \eqn{\tilde p_i=\min(1,\,m\,p_i)}. Controle le **taux d'erreur familial**
#' (probabilite d'au moins un faux positif) au niveau \eqn{\alpha}, mais est
#' **conservateur** (perte de puissance quand \eqn{m} est grand).
#'
#' @param p vecteur de p-valeurs.
#' @return p-valeurs ajustees.
#' @export
p_adjust_bonferroni <- function(p) pmin(1, p * length(p))

#' Correction de Benjamini-Hochberg (controle du FDR)
#'
#' Controle le **taux de fausses decouvertes** (proportion attendue de faux
#' positifs PARMI les rejets) au niveau \eqn{\alpha} : trie les p-valeurs, applique
#' \eqn{\tilde p_{(i)}=\min_{j\ge i} m\,p_{(j)}/j}. Bien plus **puissant** que
#' Bonferroni quand beaucoup d'hypotheses sont vraiment fausses.
#'
#' @param p vecteur de p-valeurs.
#' @return p-valeurs ajustees.
#' @export
p_adjust_bh <- function(p) {
  m <- length(p); o <- order(p); ro <- order(o)
  pmin(1, rev(cummin(rev(m / seq_len(m) * p[o]))))[ro]
}
