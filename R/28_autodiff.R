# =============================================================================
# Module 28 — Différentiation automatique en mode inverse (reverse-mode AD)
# Implémente les équations de derivations/28_autodiff.qmd. R base.
# Un mini-moteur de graphe de calcul : chaque nœud enregistre son opération sur
# une « bande » (tape) ; une seule passe arrière (règle de la chaîne) calcule
# TOUS les gradients. C'est la machinerie sous la rétropropagation (Module 12).
# =============================================================================

.ad <- new.env(parent = emptyenv()); .ad$tape <- list()
#' Réinitialise la bande d'enregistrement (avant chaque graphe)
#' @export
ad_reset <- function() .ad$tape <- list()
.ad_push <- function(nd) { .ad$tape[[length(.ad$tape) + 1L]] <- nd; nd }

#' Crée un nœud de calcul (variable enregistrée sur la bande)
#'
#' Un `adnode` est un environnement (sémantique de référence) portant sa `value`,
#' son `grad` accumulé (même forme), et une fonction `backward` qui propage le
#' gradient vers ses parents (règle de la chaîne).
#'
#' @param value scalaire, vecteur ou matrice.
#' @return objet `adnode`.
#' @export
adnode <- function(value) {
  nd <- new.env(parent = emptyenv())
  nd$value <- value; nd$grad <- value * 0
  nd$backward <- function(g) NULL
  class(nd) <- "adnode"
  .ad_push(nd)
}
as_adnode <- function(x) if (inherits(x, "adnode")) x else adnode(x)
# ramène un gradient à la forme de la cible (somme si la cible est un scalaire diffusé)
.unbroadcast <- function(g, target) if (length(target) == 1L && length(g) > 1L) sum(g) else g

#' Opérateurs élémentaires enregistrés (+, -, *, /, ^)
#' @param e1,e2 `adnode` ou numériques.
#' @return `adnode`.
#' @export
Ops.adnode <- function(e1, e2) {
  a <- as_adnode(e1); b <- as_adnode(e2); op <- .Generic
  av <- a$value; bv <- b$value
  out <- adnode(get(op)(av, bv))
  out$backward <- switch(op,
    "+" = function(g) { a$grad <- a$grad + .unbroadcast(g, av);        b$grad <- b$grad + .unbroadcast(g, bv) },
    "-" = function(g) { a$grad <- a$grad + .unbroadcast(g, av);        b$grad <- b$grad - .unbroadcast(g, bv) },
    "*" = function(g) { a$grad <- a$grad + .unbroadcast(g * bv, av);   b$grad <- b$grad + .unbroadcast(g * av, bv) },
    "/" = function(g) { a$grad <- a$grad + .unbroadcast(g / bv, av);   b$grad <- b$grad + .unbroadcast(-g * av / bv^2, bv) },
    "^" = function(g) { a$grad <- a$grad + .unbroadcast(g * bv * av^(bv - 1), av) },
    stop("operateur non supporte : ", op))
  out
}

#' Fonctions élémentaires enregistrées (exp, log, sqrt, sin, cos, tanh)
#' @param x `adnode`.
#' @param ... arguments additionnels (ignorés).
#' @return `adnode`.
#' @export
Math.adnode <- function(x, ...) {
  a <- x; av <- a$value; fn <- .Generic
  val <- get(fn)(av)
  d <- switch(fn, exp = val, log = 1 / av, sqrt = 0.5 / sqrt(av),
              sin = cos(av), cos = -sin(av), tanh = 1 - tanh(av)^2,
              stop("fonction non supportee : ", fn))
  out <- adnode(val); out$backward <- function(g) a$grad <- a$grad + g * d
  out
}

#' Somme enregistrée (réduction scalaire)
#' @param ... un unique `adnode` à sommer.
#' @param na.rm ignoré (présent pour la signature du générique).
#' @return `adnode` scalaire.
#' @export
Summary.adnode <- function(..., na.rm = FALSE) {
  if (.Generic != "sum") stop("seul `sum` est supporte")
  a <- ..1; av <- a$value
  out <- adnode(sum(av)); out$backward <- function(g) a$grad <- a$grad + g   # g scalaire diffusé
  out
}

#' Produit matriciel enregistré
#'
#' \eqn{C=AB} : \eqn{\bar A \mathrel{+}= \bar C B^\top}, \eqn{\bar B \mathrel{+}=
#' A^\top\bar C}.
#'
#' @param a,b `adnode` (ou matrices).
#' @return `adnode`.
#' @export
mm <- function(a, b) {
  a <- as_adnode(a); b <- as_adnode(b); av <- a$value; bv <- b$value
  out <- adnode(av %*% bv)
  out$backward <- function(g) { a$grad <- a$grad + g %*% t(bv); b$grad <- b$grad + t(av) %*% g }
  out
}

#' Ajoute une colonne de 1 (terme de biais) — opération enregistrée
#'
#' Pratique pour un MLP : \eqn{[\,A\;\mathbf 1\,]}. Le gradient est restreint aux
#' colonnes d'origine.
#'
#' @param a `adnode` matrice n x k.
#' @return `adnode` matrice n x (k+1).
#' @export
ad_cbind1 <- function(a) {
  a <- as_adnode(a); av <- a$value
  out <- adnode(cbind(av, 1))
  out$backward <- function(g) a$grad <- a$grad + g[, seq_len(ncol(av)), drop = FALSE]
  out
}

#' Passe arrière : calcule tous les gradients par une seule rétropropagation
#'
#' Amorce \eqn{\partial L/\partial L=1} sur le nœud de sortie, puis parcourt la
#' bande en ordre inverse (topologique) en propageant chaque gradient à ses
#' parents. Après appel, `node$grad` de chaque `adnode` contient \eqn{\partial L/
#' \partial\text{node}}.
#'
#' @param node `adnode` scalaire de sortie (la « perte »).
#' @export
backward <- function(node) {
  node$grad <- node$value * 0 + 1
  tp <- .ad$tape
  for (i in rev(seq_along(tp))) tp[[i]]$backward(tp[[i]]$grad)
  invisible(node)
}

#' Gradient d'une fonction scalaire par mode inverse
#'
#' @param f fonction `adnode -> adnode scalaire` (écrite avec les opérateurs ci-dessus).
#' @param x point (numérique) où évaluer le gradient.
#' @return le gradient \eqn{\nabla f(x)} (même forme que `x`).
#' @export
ad_grad <- function(f, x) {
  ad_reset(); xn <- adnode(x); y <- f(xn); backward(y); xn$grad
}
