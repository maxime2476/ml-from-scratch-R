# =============================================================================
# Module 15 — Interprétabilité post-hoc
# Implémente les équations de derivations/15_interpretabilite.qmd. R base.
# `predict_fn` : fonction data.frame -> vecteur numérique de prédictions.
# =============================================================================

#' Dépendance partielle (PDP, éq. 15.2)
#'
#' @param predict_fn fonction `data.frame -> numeric`.
#' @param X data.frame des données.
#' @param feature nom de la variable d'intérêt.
#' @param grid grille de valeurs (sinon régulière sur l'étendue observée).
#' @param grid_size taille de grille par défaut.
#' @return data.frame `grid`, `pdp`.
pdp <- function(predict_fn, X, feature, grid = NULL, grid_size = 25L) {
  vals <- X[[feature]]
  if (is.null(grid)) grid <- seq(min(vals), max(vals), length.out = grid_size)
  y <- vapply(grid, function(g) { Xg <- X; Xg[[feature]] <- g; mean(predict_fn(Xg)) }, numeric(1))
  data.frame(grid = grid, pdp = y)
}

#' Espérances conditionnelles individuelles (ICE, éq. 15.3)
#'
#' Le PDP est la moyenne des courbes ICE (colonne `pdp`).
#'
#' @param predict_fn fonction `data.frame -> numeric`.
#' @param X data.frame.
#' @param feature nom de la variable.
#' @param grid grille (sinon régulière).
#' @param grid_size taille par défaut.
#' @return liste : `grid`, `ice` (n x |grid|), `pdp` (moyenne des courbes).
ice <- function(predict_fn, X, feature, grid = NULL, grid_size = 25L) {
  vals <- X[[feature]]
  if (is.null(grid)) grid <- seq(min(vals), max(vals), length.out = grid_size)
  M <- vapply(grid, function(g) { Xg <- X; Xg[[feature]] <- g; predict_fn(Xg) },
              numeric(nrow(X)))
  list(grid = grid, ice = M, pdp = colMeans(M))
}

#' Valeurs de Shapley exactes par énumération (éq. 15.4, p <= 10)
#'
#' Fonction de valeur **interventionnelle** \eqn{v(S)=\mathbb E_{X_C}[f(x_S,X_C)]}
#' estimée sur `X_ref`. Complexité \eqn{O(2^p)}.
#'
#' @param predict_fn fonction `data.frame -> numeric`.
#' @param x data.frame d'une ligne (le point à expliquer).
#' @param X_ref data.frame de référence (distribution des variables absentes).
#' @return vecteur nommé des valeurs SHAP (somme = f(x) - E[f(X_ref)]).
shapley_exact <- function(predict_fn, x, X_ref) {
  features <- colnames(X_ref); p <- length(features)
  if (p > 12) stop("shapley_exact : p trop grand (utiliser shapley_permutation).")
  x <- as.data.frame(x)[features]
  bits <- 2^(0:(p - 1))
  v <- numeric(2^p)
  for (mask in 0:(2^p - 1)) {
    S <- which(bitwAnd(mask, bits) > 0)
    Z <- X_ref
    for (j in S) Z[[features[j]]] <- x[[features[j]]]
    v[mask + 1] <- mean(predict_fn(Z))
  }
  phi <- numeric(p)
  for (j in seq_len(p)) {
    bj <- bits[j]
    for (mask in 0:(2^p - 1)) {
      if (bitwAnd(mask, bj) == 0) {                       # S ne contient pas j
        s <- sum(bitwAnd(mask, bits) > 0)
        w <- factorial(s) * factorial(p - s - 1) / factorial(p)
        phi[j] <- phi[j] + w * (v[bitwOr(mask, bj) + 1] - v[mask + 1])
      }
    }
  }
  setNames(phi, features)
}

#' Valeurs de Shapley approchées par échantillonnage de permutations
#'
#' Estimateur de Štrumbelj-Kononenko : pour chaque tirage, une permutation et une
#' ligne de référence donnent la contribution marginale de chaque variable.
#'
#' @param predict_fn fonction `data.frame -> numeric`.
#' @param x data.frame d'une ligne.
#' @param X_ref data.frame de référence.
#' @param n_samples nombre de tirages.
#' @param seed graine.
#' @return vecteur nommé des valeurs SHAP estimées.
shapley_permutation <- function(predict_fn, x, X_ref, n_samples = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  features <- colnames(X_ref); p <- length(features); m <- nrow(X_ref)
  x <- as.data.frame(x)[features]
  phi <- numeric(p)
  for (s in seq_len(n_samples)) {
    perm <- sample.int(p); z <- X_ref[sample.int(m, 1), , drop = FALSE]
    xp <- z; xm <- z
    for (pos in seq_len(p)) {
      j <- perm[pos]
      xp[[features[j]]] <- x[[features[j]]]                 # ajoute j à la coalition
      phi[j] <- phi[j] + (predict_fn(xp) - predict_fn(xm))
      xm[[features[j]]] <- x[[features[j]]]                 # xm rattrape xp pour la suite
    }
  }
  setNames(phi / n_samples, features)
}

#' Importance par permutation (éq. 15.6)
#'
#' Augmentation de la perte quand on permute chaque variable (brisant son lien
#' avec y). Perte quadratique par défaut.
#'
#' @param predict_fn fonction `data.frame -> numeric`.
#' @param X data.frame.
#' @param y réponse.
#' @param loss fonction `(yhat, y) -> perte` (défaut EQM).
#' @param n_repeat répétitions de la permutation (moyennées).
#' @param seed graine.
#' @return vecteur nommé des importances.
permutation_importance <- function(predict_fn, X, y,
                                   loss = function(yh, y) mean((yh - y)^2),
                                   n_repeat = 10L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  base <- loss(predict_fn(X), y)
  imp <- vapply(colnames(X), function(feat) {
    mean(replicate(n_repeat, {
      Xp <- X; Xp[[feat]] <- sample(Xp[[feat]]); loss(predict_fn(Xp), y) - base
    }))
  }, numeric(1))
  imp
}
