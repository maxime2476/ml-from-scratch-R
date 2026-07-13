# =============================================================================
# Module 9 — Bagging et forêts aléatoires
# Implémente les équations de derivations/09_bagging_rf.qmd.
# Les arbres de base réutilisent cart_fit (Module 8), avec mtry par split pour
# la forêt aléatoire. R base + Module 8.
# =============================================================================

#' Bagging / forêt aléatoire (bootstrap + agrégation)
#'
#' Ajuste B arbres (Module 8) sur des rééchantillons bootstrap ; agrège par
#' moyenne (régression) ou vote majoritaire (classification). Avec `mtry`, chaque
#' split ne considère qu'un sous-ensemble aléatoire de variables (forêt
#' aléatoire, décorrélation des arbres). Calcule l'erreur out-of-bag (éq. 9.2).
#'
#' @param formula formule (prédicteurs numériques).
#' @param data data.frame.
#' @param method "class" ou "anova".
#' @param B nombre d'arbres.
#' @param mtry variables candidates par split ; NULL = bagging (toutes) ;
#'   sinon forêt aléatoire. Défaut : sqrt(p) (class), p/3 (anova).
#' @param max_depth,min_split,min_leaf hyperparamètres des arbres de base.
#' @param seed graine.
#' @return objet `forest` : `trees`, `oob_idx`, `oob_error`, `oob_pred`, méta.
#' @export
bagging_fit <- function(formula, data, method = c("class", "anova"), B = 100L,
                        mtry = NULL, max_depth = 30L, min_split = 5L, min_leaf = 1L,
                        seed = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  vars <- colnames(mf)[-1]
  p <- length(vars); n <- nrow(data)
  if (is.null(mtry) && method == "class") mtry <- max(1L, floor(sqrt(p)))
  if (is.null(mtry) && method == "anova") mtry <- max(1L, floor(p / 3))
  classes <- if (method == "class") levels(as.factor(y)) else NULL

  trees <- vector("list", B); oob_idx <- vector("list", B)
  # Matrice des prédictions de chaque arbre sur TOUT l'échantillon (n x B).
  pred_all <- matrix(NA_character_, n, B)
  for (b in seq_len(B)) {
    ib <- sample.int(n, n, replace = TRUE)             # échantillon bootstrap
    trees[[b]] <- cart_fit(formula, data[ib, , drop = FALSE], method = method,
                           max_depth = max_depth, min_split = min_split,
                           min_leaf = min_leaf, mtry = mtry)
    oob_idx[[b]] <- setdiff(seq_len(n), unique(ib))     # observations OOB
    pred_all[, b] <- as.character(predict_cart(trees[[b]], data))
  }

  # Agrégation OOB (éq. 9.2) : pour chaque i, moyenne/vote des arbres où i est OOB.
  in_bag <- matrix(TRUE, n, B)
  for (b in seq_len(B)) in_bag[oob_idx[[b]], b] <- FALSE   # FALSE = OOB
  oob_pred <- vector(length = n)
  for (i in seq_len(n)) {
    oobb <- which(!in_bag[i, ])
    if (length(oobb) == 0) { oob_pred[i] <- NA; next }
    pi <- pred_all[i, oobb]
    oob_pred[i] <- if (method == "class") names(which.max(table(pi)))
                   else mean(as.numeric(pi))
  }
  oob_error <- if (method == "class")
    mean(oob_pred != as.character(y), na.rm = TRUE)         # taux d'erreur
  else mean((as.numeric(oob_pred) - y)^2, na.rm = TRUE)     # EQM

  structure(list(trees = trees, oob_idx = oob_idx, oob_pred = oob_pred,
                 oob_error = oob_error, method = method, classes = classes,
                 vars = vars, B = B, mtry = mtry, n = n),
            class = "forest")
}

#' Forêt aléatoire (alias de bagging_fit avec mtry actif)
#'
#' @inheritParams bagging_fit
#' @param ... arguments supplémentaires transmis à `bagging_fit`.
#' @return objet `forest`.
#' @export
random_forest_fit <- function(formula, data, method = c("class", "anova"), B = 100L,
                              mtry = NULL, ...) {
  bagging_fit(formula, data, method = match.arg(method), B = B, mtry = mtry, ...)
}

#' Prédiction d'une forêt / d'un ensemble baggé
#'
#' Agrège les B arbres : moyenne (régression) ou vote majoritaire (classification).
#'
#' @param object objet `forest`.
#' @param newdata data.frame des prédicteurs.
#' @return vecteur de prédictions.
#' @export
predict_forest <- function(object, newdata) {
  m <- nrow(newdata); B <- object$B
  P <- matrix(NA_character_, m, B)
  for (b in seq_len(B)) P[, b] <- as.character(predict_cart(object$trees[[b]], newdata))
  if (object$method == "class") {
    out <- apply(P, 1, function(row) names(which.max(table(row))))
    factor(out, levels = object$classes)
  } else rowMeans(matrix(as.numeric(P), m, B))
}
