# =============================================================================
# Module 13 — Théorie de l'apprentissage : helpers pour les illustrations
# Implémente les objets numériques de derivations/13_slt.qmd.
# Module THÉORIQUE : ces fonctions servent à illustrer les bornes, pas à
# ajuster un modèle.
# =============================================================================

#' Borne de Hoeffding bilatérale (éq. 13.2)
#'
#' \eqn{\mathbb P(|\hat R_n - R|\ge \varepsilon)\le 2e^{-2n\varepsilon^2}} pour une
#' perte dans `[0,1]`.
#'
#' @param n taille d'échantillon.
#' @param eps écart \eqn{\varepsilon}.
#' @return la borne supérieure de probabilité.
#' @export
hoeffding_bound <- function(n, eps) pmin(1, 2 * exp(-2 * n * eps^2))

#' Complexité de Rademacher empirique d'une classe linéaire à norme bornée
#'
#' Estime \eqn{\hat{\mathfrak R}_S = \frac{B}{n}\,\mathbb E_\sigma\|\sum_i\sigma_i x_i\|}
#' (Déf. 13.10 spécialisée au cas linéaire) par tirages de signes de Rademacher.
#'
#' @param X matrice n x d des points.
#' @param B rayon de la boule \eqn{\|w\|_2\le B} (défaut 1).
#' @param n_draws nombre de tirages de sigma.
#' @param seed graine.
#' @return l'estimation de la complexité de Rademacher empirique.
#' @export
empirical_rademacher_linear <- function(X, B = 1, n_draws = 2000L, seed = NULL) {
  X <- as.matrix(X); n <- nrow(X)
  if (!is.null(seed)) set.seed(seed)
  vals <- numeric(n_draws)
  for (b in seq_len(n_draws)) {
    sigma <- sample(c(-1, 1), n, replace = TRUE)
    vals[b] <- sqrt(sum((crossprod(X, sigma))^2))    # ||sum sigma_i x_i||
  }
  B * mean(vals) / n
}

#' Borne théorique de Rademacher pour une classe linéaire (éq. 13.5)
#'
#' \eqn{\hat{\mathfrak R}_S \le B\rho/\sqrt n}, \eqn{\rho=\max_i\|x_i\|}.
#'
#' @param X matrice n x d des points.
#' @param B rayon de la boule (défaut 1).
#' @return la borne \eqn{B\rho/\sqrt n}.
#' @export
rademacher_linear_bound <- function(X, B = 1) {
  X <- as.matrix(X)
  rho <- max(sqrt(rowSums(X^2)))
  B * rho / sqrt(nrow(X))
}

#' Distance de l'origine à l'enveloppe convexe de colonnes (QP exact)
#'
#' Résout \eqn{\min_{\lambda\ge 0,\ \mathbf 1^\top\lambda=1}\|V\lambda\|} par
#' programmation quadratique. Sert au test de séparabilité (théorème de Gordan :
#' un étiquetage est linéairement séparable ssi 0 n'est PAS dans l'enveloppe
#' convexe des \eqn{y_i\,\tilde x_i}).
#'
#' @param V matrice (colonnes = points).
#' @return la distance minimale de l'origine à l'enveloppe convexe.
#' @keywords internal
.min_dist_hull <- function(V) {
  V <- as.matrix(V); m <- ncol(V)
  Dmat <- crossprod(V) + diag(1e-10, m)          # V'V (régularisée pour la SDP)
  dvec <- rep(0, m)
  Amat <- cbind(rep(1, m), diag(m))               # Sum lambda = 1 ; lambda >= 0
  bvec <- c(1, rep(0, m))
  sol <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  lam <- pmax(sol$solution, 0); lam <- lam / sum(lam)
  sqrt(sum((V %*% lam)^2))
}

#' Séparabilité linéaire d'un étiquetage (théorème de Gordan)
#'
#' @param X matrice n x d des points.
#' @param y étiquettes dans \{-1, +1\} (longueur n).
#' @param tol seuil de distance pour déclarer 0 hors de l'enveloppe.
#' @return TRUE si \{(x_i, y_i)\} est linéairement séparable (avec biais).
#' @export
is_separable <- function(X, y, tol = 1e-6) {
  Xa <- cbind(as.matrix(X), 1)                    # augmentation (biais)
  V <- t(Xa * y)                                   # colonnes : y_i * aug(x_i)
  .min_dist_hull(V) > tol                          # 0 hors enveloppe -> séparable
}

#' La classe des hyperplans pulvérise-t-elle un ensemble de points ? (Déf. 13.6)
#'
#' Teste si les \eqn{2^m} étiquetages d'un ensemble de m points sont TOUS
#' linéairement séparables (donc réalisables par un demi-espace).
#'
#' @param X matrice m x d des points.
#' @return TRUE si l'ensemble est pulvérisé par les hyperplans de R^d.
#' @export
shatters_hyperplane <- function(X) {
  X <- as.matrix(X); m <- nrow(X)
  labelings <- expand.grid(rep(list(c(-1, 1)), m))
  for (r in seq_len(nrow(labelings))) {
    y <- as.numeric(labelings[r, ])
    if (length(unique(y)) == 1) next               # étiquetage constant : trivial
    if (!is_separable(X, y)) return(FALSE)
  }
  TRUE
}
