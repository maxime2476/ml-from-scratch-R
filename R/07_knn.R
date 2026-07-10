# =============================================================================
# Module 7 — K plus proches voisins et fléau de la dimension
# Implémente les équations de derivations/07_knn.qmd. R base uniquement.
# =============================================================================

# Distances euclidiennes de chaque ligne de Xte à toutes les lignes de Xtr.
.dist_to_train <- function(Xtr, xrow) sqrt(colSums((t(Xtr) - xrow)^2))

#' Régression KNN (éq. 7.1)
#'
#' Prédit en chaque point de test la moyenne des réponses de ses k plus proches
#' voisins d'apprentissage.
#'
#' @param X_train matrice n x p d'apprentissage.
#' @param y_train réponses d'apprentissage (numériques).
#' @param X_test matrice m x p de test.
#' @param k nombre de voisins.
#' @return vecteur des m prédictions.
knn_regression <- function(X_train, y_train, X_test, k) {
  Xtr <- as.matrix(X_train); Xte <- as.matrix(X_test)
  if (k > nrow(Xtr)) stop("k > nombre d'observations d'apprentissage.")
  vapply(seq_len(nrow(Xte)), function(i) {
    d <- .dist_to_train(Xtr, Xte[i, ])
    mean(y_train[order(d)[seq_len(k)]])
  }, numeric(1))
}

#' Classification KNN par vote majoritaire
#'
#' Vote majoritaire des k plus proches voisins. Les égalités de vote sont
#' rompues par la première classe (ordre de `levels`) — pour une comparaison
#' déterministe avec `class::knn`, utiliser k impair et 2 classes (pas d'égalité).
#'
#' @param X_train matrice n x p d'apprentissage.
#' @param y_train étiquettes d'apprentissage (facteur ou vecteur).
#' @param X_test matrice m x p de test.
#' @param k nombre de voisins.
#' @return facteur des m classes prédites (mêmes niveaux que y_train).
knn_classify <- function(X_train, y_train, X_test, k) {
  Xtr <- as.matrix(X_train); Xte <- as.matrix(X_test)
  y <- as.factor(y_train); lev <- levels(y)
  if (k > nrow(Xtr)) stop("k > nombre d'observations d'apprentissage.")
  out <- vapply(seq_len(nrow(Xte)), function(i) {
    d <- .dist_to_train(Xtr, Xte[i, ])
    votes <- table(y[order(d)[seq_len(k)]])
    names(votes)[which.max(votes)]              # which.max : 1ère classe si égalité
  }, character(1))
  factor(out, levels = lev)
}

#' Longueur d'arête d'un voisinage capturant une fraction r (éq. 7.3)
#'
#' \eqn{e_p(r) = r^{1/p}} dans l'hypercube unité.
#'
#' @param r fraction des données à capturer (dans (0,1]).
#' @param p dimension.
#' @return la longueur d'arête nécessaire.
edge_length <- function(r, p) r^(1 / p)

#' Concentration des distances (éq. 7.4-7.5)
#'
#' Pour n points i.i.d. en dimension p, calcule, pour chaque point de requête, le
#' rapport de contraste \eqn{(D_{\max}-D_{\min})/D_{\min}} entre distance la plus
#' grande et la plus petite aux autres points, et le coefficient de variation
#' des distances au carré (théorique \eqn{\propto 1/\sqrt p}).
#'
#' @param n nombre de points.
#' @param p dimension.
#' @param gen générateur de points (défaut : uniforme sur [0,1]^p).
#' @param seed graine.
#' @return liste : `contrast` (moyenne de (Dmax-Dmin)/Dmin), `cv_d2` (coef. de
#'   variation des distances au carré).
distance_concentration <- function(n, p, gen = function(n, p) matrix(runif(n * p), n, p),
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- gen(n, p)
  D <- as.matrix(dist(X))                       # distances par paires
  diag(D) <- NA
  dmin <- apply(D, 1, min, na.rm = TRUE)
  dmax <- apply(D, 1, max, na.rm = TRUE)
  d2 <- D[upper.tri(D)]^2
  list(contrast = mean((dmax - dmin) / dmin),
       cv_d2 = sd(d2) / mean(d2))
}
