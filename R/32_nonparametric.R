# =============================================================================
# Module 32 — Regression non parametrique et discontinuite de regression (RDD)
# Implemente les equations de derivations/32_nonparametric.qmd. R base.
# Au lieu d'imposer une forme (lineaire, ...), on estime m(x)=E[Y|X=x] LOCALEMENT
# par ponderation a noyau. Application : la discontinuite de regression (RDD),
# design causal quasi-experimental fonde sur la regression locale.
# =============================================================================

.gauss <- function(u) dnorm(u)
.tri <- function(u) pmax(0, 1 - abs(u))                    # noyau triangulaire

#' Estimation de densite par noyau (KDE)
#'
#' \eqn{\hat f(x)=\frac1{nh}\sum_i K\!\bigl(\frac{x-x_i}{h}\bigr)}, noyau gaussien.
#'
#' @param x donnees ; @param grid points d'evaluation ; @param bw fenetre h.
#' @return vecteur des densites estimees sur `grid`.
#' @export
kde <- function(x, grid, bw) {
  vapply(grid, function(g) mean(.gauss((g - x) / bw)) / bw, numeric(1))
}

#' Regression de Nadaraya-Watson (local constant)
#'
#' \eqn{\hat m(x_0)=\sum_i K_h(x_0-x_i)\,y_i/\sum_i K_h(x_0-x_i)} : moyenne locale
#' ponderee par le noyau. Biais d'ordre \eqn{h^2}, mais biais de bord marque.
#'
#' @param x,y donnees ; @param x0 points d'evaluation ; @param bw fenetre.
#' @return vecteur des valeurs ajustees en `x0`.
#' @export
nadaraya_watson <- function(x, y, x0, bw) {
  vapply(x0, function(g) { w <- .gauss((g - x) / bw); sum(w * y) / sum(w) }, numeric(1))
}

#' Regression locale lineaire
#'
#' Ajuste en chaque \eqn{x_0} une droite ponderee par le noyau (WLS locale) et
#' renvoie l'ordonnee a l'origine. Corrige le **biais de bord** du Nadaraya-Watson
#' (biais d'ordre \eqn{h^2} uniforme jusqu'aux bords).
#'
#' @param x,y donnees ; @param x0 points d'evaluation ; @param bw fenetre ;
#' @param kernel "gauss" ou "tri" (triangulaire).
#' @return vecteur des valeurs ajustees en `x0`.
#' @export
local_linear <- function(x, y, x0, bw, kernel = c("gauss", "tri")) {
  K <- if (match.arg(kernel) == "gauss") .gauss else .tri
  vapply(x0, function(g) {
    w <- K((g - x) / bw); if (sum(w > 0) < 2) return(NA_real_)
    X <- cbind(1, x - g); b <- solve(crossprod(X * w, X), crossprod(X * w, y))
    b[1]
  }, numeric(1))
}

#' Selection de fenetre par validation croisee leave-one-out (Nadaraya-Watson)
#'
#' Minimise \eqn{\sum_i (y_i-\hat m_{-i}(x_i))^2}, ou \eqn{\hat m_{-i}} exclut
#' l'observation \eqn{i} (formule du residu : retirer le poids diagonal).
#'
#' @param x,y donnees ; @param bws fenetres candidates.
#' @return liste : `bw` (optimale), `cv` (erreurs CV par fenetre).
#' @export
bw_loocv <- function(x, y, bws) {
  cv <- vapply(bws, function(h) {
    err <- vapply(seq_along(x), function(i) {
      w <- .gauss((x[i] - x[-i]) / h); (y[i] - sum(w * y[-i]) / sum(w))^2
    }, numeric(1))
    mean(err)
  }, numeric(1))
  list(bw = bws[which.min(cv)], cv = cv)
}

#' Discontinuite de regression (RDD, effet local au seuil)
#'
#' Ajuste une regression **locale lineaire** de chaque cote d'un seuil \eqn{c} et
#' estime le **saut** \eqn{\tau=\lim_{x\downarrow c}m(x)-\lim_{x\uparrow c}m(x)},
#' l'effet causal local (design quasi-experimental de Thistlethwaite-Campbell).
#'
#' @param y resultat ; @param x variable de forcage (running variable).
#' @param cutoff seuil \eqn{c} ; @param bw fenetre (des deux cotes).
#' @param kernel "tri" (triangulaire, defaut) ou "gauss".
#' @return liste : `tau` (saut), `mu_left`, `mu_right`, `n_left`, `n_right`.
#' @export
rdd <- function(y, x, cutoff = 0, bw, kernel = c("tri", "gauss")) {
  K <- if (match.arg(kernel) == "tri") .tri else .gauss
  side <- function(idx) {
    w <- K((x[idx] - cutoff) / bw); X <- cbind(1, x[idx] - cutoff)
    as.numeric(solve(crossprod(X * w, X), crossprod(X * w, y[idx]))[1])
  }
  L <- which(x < cutoff & x >= cutoff - bw); R <- which(x >= cutoff & x <= cutoff + bw)
  muL <- side(L); muR <- side(R)
  list(tau = muR - muL, mu_left = muL, mu_right = muR, n_left = length(L), n_right = length(R))
}
