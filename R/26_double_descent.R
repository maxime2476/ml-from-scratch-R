# =============================================================================
# Module 26 — Double descente (régime sur-paramétré)
# Implémente les équations de derivations/26_double_descent.qmd.
# Réutilise solve_ls_svd (M0, interpolateur de norme minimale) et ridge_fit (M4).
# =============================================================================

#' Caractéristiques de Fourier aléatoires (approximation du noyau gaussien)
#'
#' \eqn{\varphi(x)=\sqrt{2/D}\,\cos(Wx+b)}, avec \eqn{W\sim\mathcal N(0,\gamma I)}
#' et \eqn{b\sim\mathcal U(0,2\pi)}. La dimension \eqn{D} des caractéristiques est
#' le **paramètre de complexité** dont on fera varier la valeur autour de \eqn{n}.
#'
#' @param X matrice n x p.
#' @param D nombre de caractéristiques.
#' @param gamma échelle (largeur inverse du noyau).
#' @param seed graine (les poids aléatoires doivent être FIXES entre train/test).
#' @return matrice n x D.
#' @export
random_features <- function(X, D, gamma = 1, seed = 1) {
  X <- as.matrix(X); p <- ncol(X)
  old <- .Random.seed; on.exit({ .Random.seed <<- old }, add = TRUE); set.seed(seed)
  W <- matrix(rnorm(p * D, sd = sqrt(gamma)), p, D)
  b <- runif(D, 0, 2 * pi)
  sqrt(2 / D) * cos(sweep(X %*% W, 2, b, "+"))
}

#' Ajuste un modèle à caractéristiques aléatoires (interpolant ou régularisé)
#'
#' `lambda = 0` : moindres carrés de **norme minimale** (`solve_ls_svd`, M0) —
#' interpole les données dès que \eqn{D \ge n}. `lambda > 0` : ridge (M4) sur les
#' caractéristiques. Renvoie un prédicteur.
#'
#' @param X,y données d'apprentissage.
#' @param D dimension des caractéristiques ; @param gamma,seed cf. `random_features`.
#' @param lambda pénalité ridge (0 = interpolation de norme minimale).
#' @return fonction `newX -> prédictions`.
#' @export
fit_rff <- function(X, y, D, gamma = 1, seed = 1, lambda = 0) {
  Phi <- random_features(X, D, gamma, seed)
  if (lambda == 0) {
    beta <- solve_ls_svd(Phi, y)$coefficients            # norme minimale (M0)
    function(newX) as.numeric(random_features(newX, D, gamma, seed) %*% beta)
  } else {
    fit <- ridge_fit(Phi, y, lambda = lambda, standardize = FALSE, intercept = TRUE)  # M4
    function(newX) {
      P <- random_features(newX, D, gamma, seed)
      fit$intercept + as.numeric(P %*% fit$beta)          # beta = pentes (sans intercept)
    }
  }
}

#' Courbe de risque en fonction de la complexité (double descente)
#'
#' Balaie la dimension \eqn{D} des caractéristiques et renvoie l'erreur
#' quadratique d'**apprentissage** et de **test** pour chaque \eqn{D}. Avec
#' `lambda = 0`, le test présente un **pic à \eqn{D\approx n}** (seuil
#' d'interpolation) puis **redescend** — la double descente.
#'
#' @param Xtr,ytr apprentissage ; @param Xte,yte test.
#' @param Ds vecteur des dimensions à essayer.
#' @param gamma,seed cf. `random_features` ; @param lambda pénalité ridge.
#' @return data.frame : `D`, `train_mse`, `test_mse`.
#' @export
double_descent_curve <- function(Xtr, ytr, Xte, yte, Ds, gamma = 1, seed = 1, lambda = 0) {
  do.call(rbind, lapply(Ds, function(D) {
    pred <- fit_rff(Xtr, ytr, D, gamma, seed, lambda)
    data.frame(D = D,
               train_mse = mean((ytr - pred(Xtr))^2),
               test_mse  = mean((yte - pred(Xte))^2))
  }))
}
