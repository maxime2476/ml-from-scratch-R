# =============================================================================
# Module 0 — Algèbre linéaire numérique et optimiseurs génériques
# Implémente les équations de derivations/00_linalg.qmd.
#
# Rappel du parti pris (cf. README, "Auto-suffisance par module") : les
# optimiseurs optim_*() de ce fichier sont la VERSION CANONIQUE de référence
# (et le banc d'essai des Monte Carlo). Les modules suivants réimplémentent
# localement leur propre optimiseur spécialisé plutôt que d'importer ceux-ci.
# =============================================================================

# ---- Solveurs triangulaires (briques de base, from scratch) -----------------

#' Résolution d'un système triangulaire inférieur par descente
#'
#' Résout \eqn{L x = b} avec \eqn{L} triangulaire inférieure. Brique de la
#' résolution des équations normales par Cholesky, éq. (0.10) de
#' derivations/00_linalg.qmd.
#'
#' @param L matrice triangulaire inférieure (n x n), diagonale non nulle.
#' @param b vecteur second membre (longueur n).
#' @return le vecteur solution x.
#' @export
forward_substitution <- function(L, b) {
  n <- length(b)
  if (nrow(L) != n || ncol(L) != n) stop("Dimensions L/b incompatibles.")
  x <- numeric(n)
  for (i in seq_len(n)) {
    prev <- seq_len(i - 1L)
    x[i] <- (b[i] - sum(L[i, prev] * x[prev])) / L[i, i]
  }
  x
}

#' Résolution d'un système triangulaire supérieur par remontée
#'
#' Résout \eqn{U x = b} avec \eqn{U} triangulaire supérieure. Brique de la
#' résolution MCO par QR (remontée sur \eqn{R_1}, éq. 0.8) et par Cholesky
#' (éq. 0.10).
#'
#' @param U matrice triangulaire supérieure (n x n), diagonale non nulle.
#' @param b vecteur second membre (longueur n).
#' @return le vecteur solution x.
#' @export
back_substitution <- function(U, b) {
  n <- length(b)
  if (nrow(U) != n || ncol(U) != n) stop("Dimensions U/b incompatibles.")
  x <- numeric(n)
  for (i in rev(seq_len(n))) {
    nxt <- if (i < n) (i + 1L):n else integer(0)
    x[i] <- (b[i] - sum(U[i, nxt] * x[nxt])) / U[i, i]
  }
  x
}

# ---- QR par réflexions de Householder ---------------------------------------

#' Vecteur de Householder (choix de signe stable)
#'
#' Renvoie le vecteur \eqn{v} définissant la réflexion \eqn{H(v)=I-2vv^T/v^Tv}
#' telle que \eqn{H(v)x = \alpha e_1}. Implémente l'éq. (0.5) : le signe
#' \eqn{v_1 = x_1 + \mathrm{sign}(x_1)\|x\|} évite l'annulation catastrophique.
#'
#' @param x vecteur non nul.
#' @return le vecteur v (même longueur que x) ; le vecteur nul si x est nul.
#' @export
householder_vector <- function(x) {
  nx <- sqrt(sum(x^2))
  if (nx == 0) return(rep(0, length(x)))
  s <- if (x[1] >= 0) 1 else -1
  v <- x
  v[1] <- v[1] + s * nx
  v
}

#' Décomposition QR par réflexions de Householder
#'
#' Factorise \eqn{X = QR} avec \eqn{Q} orthogonale (n x n) et \eqn{R}
#' triangulaire supérieure (n x p). Implémente les éq. (0.4)-(0.6) de
#' derivations/00_linalg.qmd. \eqn{Q} est accumulée explicitement
#' (\eqn{Q = H_1\cdots H_p}) à des fins pédagogiques.
#'
#' @param X matrice n x p, \eqn{n \ge p}.
#' @return liste `Q` (n x n orthogonale), `R` (n x p triangulaire supérieure).
#' @export
qr_householder <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (n < p) stop("qr_householder attend n >= p.")
  R <- X
  Q <- diag(n)
  for (j in seq_len(min(p, n - 1L))) {
    idx <- j:n
    v <- householder_vector(R[idx, j])
    vtv <- sum(v * v)
    if (vtv == 0) next
    # Applique H_j à gauche sur R[idx, j:p] : R <- R - (2/v'v) v (v' R)
    R[idx, j:p] <- R[idx, j:p] - (2 / vtv) * (v %*% (t(v) %*% R[idx, j:p]))
    # Accumule Q à droite : Q[, idx] <- Q[, idx] - (2/v'v) (Q[, idx] v) v'
    Q[, idx] <- Q[, idx] - (2 / vtv) * ((Q[, idx] %*% v) %*% t(v))
  }
  # Nettoie la poussière numérique sous la diagonale de R.
  for (j in seq_len(p)) if (j < n) R[(j + 1L):n, j] <- 0
  list(Q = Q, R = R)
}

#' Moindres carrés par QR de Householder
#'
#' Résout \eqn{\min_\beta \|X\beta - y\|^2} via la factorisation QR, en
#' exploitant la préservation de norme par \eqn{Q} : éq. (0.7)-(0.8). La somme
#' des carrés des résidus est obtenue « gratuitement » comme \eqn{\|c_2\|^2}.
#' Voie recommandée quand X est bien conditionnée (préserve \eqn{\kappa_2(X)}).
#'
#' @param X matrice de design n x p, de plein rang colonne.
#' @param y vecteur réponse (longueur n).
#' @return liste : `coefficients`, `fitted`, `residuals`, `rss`, `R` (\eqn{R_1}),
#'   `Qty` (\eqn{Q^T y}).
#' @export
solve_ls_qr <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  if (length(y) != n) stop("length(y) doit valoir nrow(X).")
  dec <- qr_householder(X)
  Qty <- as.numeric(crossprod(dec$Q, y))          # Q^T y
  R1  <- dec$R[seq_len(p), , drop = FALSE]          # bloc p x p
  c1  <- Qty[seq_len(p)]
  c2  <- if (n > p) Qty[(p + 1L):n] else numeric(0)
  dR <- abs(diag(R1))
  if (any(dR < .Machine$double.eps * max(dR) * max(n, p)))
    stop("X colinéaire (R1 singulière) : utiliser solve_ls_svd().")
  beta <- back_substitution(R1, c1)
  fitted <- as.numeric(X %*% beta)
  list(coefficients = beta, fitted = fitted, residuals = y - fitted,
       rss = sum(c2^2), R = R1, Qty = Qty)
}

# ---- Factorisation de Cholesky ----------------------------------------------

#' Factorisation de Cholesky (algorithme de Crout)
#'
#' Factorise une matrice SPD \eqn{A = L L^T}, \eqn{L} triangulaire inférieure à
#' diagonale positive. Implémente les formules (0.9) de derivations/00_linalg.qmd.
#'
#' @param A matrice carrée symétrique définie positive.
#' @return le facteur `L` (triangulaire inférieur).
#' @export
chol_crout <- function(A) {
  A <- as.matrix(A)
  p <- nrow(A)
  if (ncol(A) != p) stop("A doit être carrée.")
  L <- matrix(0, p, p)
  for (j in seq_len(p)) {
    prev <- seq_len(j - 1L)
    s <- A[j, j] - sum(L[j, prev]^2)
    if (s <= 0) stop("A non définie positive (pivot <= 0).")
    L[j, j] <- sqrt(s)
    if (j < p) for (i in (j + 1L):p) {
      L[i, j] <- (A[i, j] - sum(L[i, prev] * L[j, prev])) / L[j, j]
    }
  }
  L
}

#' Moindres carrés par équations normales et Cholesky
#'
#' Résout \eqn{X^T X \beta = X^T y} en factorisant \eqn{X^T X = L L^T} puis par
#' descente/remontée : éq. (0.10). Rapide mais subit \eqn{\kappa_2(X)^2}
#' (cf. Prop. 0.1) — à réserver aux problèmes bien conditionnés.
#'
#' @param X matrice de design n x p, de plein rang colonne.
#' @param y vecteur réponse (longueur n).
#' @return liste : `coefficients`, `fitted`, `residuals`, `rss`, `L`.
#' @export
solve_ls_chol <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y)
  if (length(y) != nrow(X)) stop("length(y) doit valoir nrow(X).")
  XtX <- crossprod(X)                 # X^T X
  Xty <- as.numeric(crossprod(X, y))  # X^T y
  L <- chol_crout(XtX)
  z    <- forward_substitution(L, Xty)     # L z = X^T y
  beta <- back_substitution(t(L), z)        # L^T beta = z
  fitted <- as.numeric(X %*% beta)
  list(coefficients = beta, fitted = fitted, residuals = y - fitted,
       rss = sum((y - fitted)^2), L = L)
}

# ---- SVD : rang, conditionnement, pseudo-inverse ----------------------------

#' Outils SVD : rang numérique, conditionnement, pseudo-inverse
#'
#' À partir de \eqn{X = U\Sigma V^T} (éq. 0.11), calcule le rang numérique, le
#' conditionnement \eqn{\kappa_2 = \sigma_{\max}/\sigma_{\min}} (éq. 0.3) et la
#' pseudo-inverse de Moore-Penrose \eqn{X^+ = V\Sigma^+ U^T} (éq. 0.12). La SVD
#' elle-même est déléguée à `svd()` (existence admise, Th. 0.4).
#'
#' @param X matrice n x p.
#' @param tol seuil de rang ; par défaut \eqn{\max(n,p)\,\varepsilon\,\sigma_1}.
#' @return liste : `d` (valeurs singulières), `rank`, `kappa`, `pinv`
#'   (pseudo-inverse p x n), `u`, `v`, `tol`.
#' @export
svd_tools <- function(X, tol = NULL) {
  X <- as.matrix(X)
  sv <- svd(X)
  d <- sv$d
  if (is.null(tol)) tol <- max(dim(X)) * .Machine$double.eps * d[1]
  pos <- d > tol
  r <- sum(pos)
  dinv <- ifelse(pos, 1 / d, 0)
  pinv <- sv$v %*% (dinv * t(sv$u))   # V Sigma+ U^T
  kappa <- if (d[length(d)] > tol) d[1] / d[length(d)] else Inf
  list(d = d, rank = r, kappa = kappa, pinv = pinv,
       u = sv$u, v = sv$v, tol = tol)
}

#' Moindres carrés de norme minimale par SVD
#'
#' Renvoie \eqn{\hat\beta_{\min} = X^+ y} (éq. 0.13), qui minimise
#' \eqn{\|X\beta - y\|} et, en cas de rang déficient, a la plus petite norme
#' parmi tous les minimiseurs (Prop. 0.5). Solution propre à la colinéarité
#' parfaite et au cas \eqn{p > n}.
#'
#' @param X matrice de design n x p (rang quelconque).
#' @param y vecteur réponse (longueur n).
#' @return liste : `coefficients`, `fitted`, `residuals`, `rss`, `rank`, `kappa`.
#' @export
solve_ls_svd <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y)
  if (length(y) != nrow(X)) stop("length(y) doit valoir nrow(X).")
  s <- svd_tools(X)
  beta <- as.numeric(s$pinv %*% y)
  fitted <- as.numeric(X %*% beta)
  list(coefficients = beta, fitted = fitted, residuals = y - fitted,
       rss = sum((y - fitted)^2), rank = s$rank, kappa = s$kappa)
}

# ---- Optimiseurs génériques (versions canoniques de référence) --------------

#' Descente de gradient à pas constant
#'
#' Itère \eqn{x_{k+1} = x_k - t\,\nabla f(x_k)} (éq. 0.15). Sous \eqn{f} convexe
#' L-lisse et pas \eqn{t = 1/L}, vitesse \eqn{O(1/k)} (Th. 0.6).
#'
#' @param grad fonction gradient : `grad(x)` renvoie \eqn{\nabla f(x)}.
#' @param x0 point initial.
#' @param step pas constant \eqn{t} (idéalement \eqn{1/L}).
#' @param max_iter nombre maximal d'itérations.
#' @param tol seuil d'arrêt sur \eqn{\|x_{k+1}-x_k\|}.
#' @param f (optionnel) fonction objectif, pour renvoyer la valeur finale.
#' @return liste : `par`, `iter`, `grad_norm`, `value`.
#' @export
optim_gd <- function(grad, x0, step, max_iter = 1e4L, tol = 1e-8, f = NULL) {
  x <- x0
  k <- 0L
  for (k in seq_len(max_iter)) {
    g <- grad(x)
    x_new <- x - step * g
    if (sqrt(sum((x_new - x)^2)) < tol) { x <- x_new; break }
    x <- x_new
  }
  list(par = x, iter = k, grad_norm = sqrt(sum(grad(x)^2)),
       value = if (!is.null(f)) f(x) else NA_real_)
}

#' Newton-Raphson (minimisation)
#'
#' Itère \eqn{x_{k+1} = x_k - [\nabla^2 f(x_k)]^{-1}\nabla f(x_k)} (éq. 0.17).
#' Convergence quadratique locale (Th. 0.7). Le pas de Newton est obtenu en
#' résolvant \eqn{\nabla^2 f\,\delta = \nabla f}.
#'
#' @param grad fonction gradient `grad(x)`.
#' @param hess fonction hessienne `hess(x)` (matrice SPD localement).
#' @param x0 point initial.
#' @param max_iter nombre maximal d'itérations.
#' @param tol seuil d'arrêt sur \eqn{\|x_{k+1}-x_k\|}.
#' @param f (optionnel) fonction objectif.
#' @return liste : `par`, `iter`, `grad_norm`, `value`.
#' @export
optim_newton <- function(grad, hess, x0, max_iter = 100L, tol = 1e-10, f = NULL) {
  x <- x0
  k <- 0L
  for (k in seq_len(max_iter)) {
    g <- grad(x); H <- hess(x)
    delta <- solve(H, g)
    x_new <- x - delta
    if (sqrt(sum((x_new - x)^2)) < tol) { x <- x_new; break }
    x <- x_new
  }
  list(par = x, iter = k, grad_norm = sqrt(sum(grad(x)^2)),
       value = if (!is.null(f)) f(x) else NA_real_)
}

#' Coordinate descent générique
#'
#' Minimise cycliquement selon chaque coordonnée (éq. 0.19), en supposant fourni
#' le minimiseur 1-D exact `argmin_coord(x, j)` (structure « sous-problème 1-D
#' fermé », qui deviendra le soft-thresholding au Module 4).
#'
#' @param argmin_coord fonction : `argmin_coord(x, j)` renvoie la valeur optimale
#'   de la coordonnée j, les autres coordonnées de x étant fixées.
#' @param x0 point initial.
#' @param max_sweep nombre maximal de balayages complets.
#' @param tol seuil d'arrêt sur la variation d'un balayage.
#' @param f (optionnel) fonction objectif.
#' @return liste : `par`, `sweeps`, `value`.
#' @export
optim_cd <- function(argmin_coord, x0, max_sweep = 1000L, tol = 1e-9, f = NULL) {
  x <- x0; d <- length(x)
  s <- 0L
  for (s in seq_len(max_sweep)) {
    x_old <- x
    for (j in seq_len(d)) x[j] <- argmin_coord(x, j)
    if (sqrt(sum((x - x_old)^2)) < tol) break
  }
  list(par = x, sweeps = s, value = if (!is.null(f)) f(x) else NA_real_)
}

#' Gradient stochastique (SGD) par mini-lots
#'
#' Minimise une perte-somme \eqn{f(x)=\frac1n\sum_i f_i(x)} par l'itération
#' \eqn{x_{k+1}=x_k - t_k g_k} avec \eqn{g_k} gradient sur un mini-lot
#' (éq. 0.20). Pas constant (`step`) ou décroissant (`step_fun`, conditions de
#' Robbins-Monro 0.21).
#'
#' @param grad_i fonction : `grad_i(x, idx)` renvoie le gradient moyen sur les
#'   observations d'indices `idx`.
#' @param x0 point initial.
#' @param n nombre total d'observations.
#' @param batch taille de mini-lot.
#' @param step pas constant (ignoré si `step_fun` fourni).
#' @param step_fun (optionnel) fonction du compteur d'updates t -> pas.
#' @param epochs nombre de passages sur l'échantillon.
#' @param seed (optionnel) graine pour la permutation des indices.
#' @return liste : `par`, `updates`.
#' @export
optim_sgd <- function(grad_i, x0, n, batch = 1L, step = NULL, step_fun = NULL,
                      epochs = 50L, seed = NULL) {
  if (is.null(step) && is.null(step_fun)) stop("Fournir step ou step_fun.")
  if (!is.null(seed)) set.seed(seed)
  x <- x0; t <- 0L
  for (e in seq_len(epochs)) {
    idx <- sample.int(n)
    for (start in seq(1L, n, by = batch)) {
      t <- t + 1L
      bi <- idx[start:min(start + batch - 1L, n)]
      g <- grad_i(x, bi)
      st <- if (!is.null(step_fun)) step_fun(t) else step
      x <- x - st * g
    }
  }
  list(par = x, updates = t)
}

#' Gradient accéléré de Nesterov
#'
#' Ajoute un terme d'inertie (« momentum ») à la descente de gradient : pour f
#' convexe L-lisse et pas \eqn{t=1/L}, la vitesse passe de \eqn{O(1/k)}
#' (gradient) à \eqn{O(1/k^2)} — l'accélération optimale au premier ordre
#' (Nesterov 1983).
#'
#' @param grad fonction gradient `grad(x)`.
#' @param x0 point initial.
#' @param step pas \eqn{t=1/L}.
#' @param max_iter itérations maximales.
#' @param tol seuil d'arrêt sur \eqn{\|x_{k+1}-x_k\|}.
#' @param f (optionnel) fonction objectif.
#' @return liste : `par`, `iter`, `grad_norm`, `value`.
#' @export
optim_nesterov <- function(grad, x0, step, max_iter = 1e4L, tol = 1e-8, f = NULL) {
  x <- x0; y <- x0; lam <- 1; k <- 0L
  for (k in seq_len(max_iter)) {
    g <- grad(y)
    x_new <- y - step * g
    # restart adaptatif (O'Donoghue-Candès) : si l'on avance CONTRE le gradient,
    # réinitialiser l'inertie -> évite l'oscillation en régime fortement convexe.
    if (sum(g * (x_new - x)) > 0) lam <- 1
    lam_new <- (1 + sqrt(1 + 4 * lam^2)) / 2
    y <- x_new + ((lam - 1) / lam_new) * (x_new - x)      # extrapolation d'inertie
    if (sqrt(sum((x_new - x)^2)) < tol) { x <- x_new; break }
    x <- x_new; lam <- lam_new
  }
  list(par = x, iter = k, grad_norm = sqrt(sum(grad(x)^2)),
       value = if (!is.null(f)) f(x) else NA_real_)
}

# Recherche linéaire de Wolfe (faible) par bissection : renvoie un pas alpha
# satisfaisant Armijo (c1) et la condition de courbure (c2).
.wolfe_ls <- function(f, grad, x, d, g0, fx, c1 = 1e-4, c2 = 0.4, max_ls = 50L) {
  slope0 <- sum(g0 * d); alpha <- 1; alo <- 0; ahi <- Inf
  for (it in seq_len(max_ls)) {
    if (f(x + alpha * d) > fx + c1 * alpha * slope0) {          # Armijo violé
      ahi <- alpha; alpha <- (alo + ahi) / 2
    } else if (sum(grad(x + alpha * d) * d) < c2 * slope0) {    # courbure violée
      alo <- alpha; alpha <- if (is.finite(ahi)) (alo + ahi) / 2 else 2 * alpha
    } else return(alpha)
  }
  alpha
}

#' L-BFGS (quasi-Newton à mémoire limitée)
#'
#' Approxime l'inverse de la hessienne à partir des \eqn{m} derniers couples
#' \eqn{(s_k,y_k)} par la **récursion à deux boucles**, sans stocker de matrice
#' \eqn{d\times d}. Recherche linéaire d'Armijo (rétrogression) si `f` est fourni.
#' Convergence super-linéaire, coût par itération \eqn{O(md)}.
#'
#' @param grad fonction gradient `grad(x)`.
#' @param x0 point initial.
#' @param f fonction objectif (pour la recherche linéaire ; recommandé).
#' @param m taille de mémoire (défaut 10).
#' @param max_iter itérations maximales.
#' @param tol seuil sur la norme du gradient.
#' @return liste : `par`, `iter`, `grad_norm`, `value`.
#' @export
optim_lbfgs <- function(grad, x0, f = NULL, m = 10L, max_iter = 200L, tol = 1e-8) {
  x <- x0; g <- grad(x); S <- list(); Y <- list(); k <- 0L
  for (k in seq_len(max_iter)) {
    if (sqrt(sum(g^2)) < tol) break
    # --- récursion à deux boucles : d = -H_k g ---
    q <- g; ns <- length(S); alphas <- numeric(ns); rho <- numeric(ns)
    for (i in rev(seq_len(ns))) {
      rho[i] <- 1 / sum(Y[[i]] * S[[i]]); alphas[i] <- rho[i] * sum(S[[i]] * q)
      q <- q - alphas[i] * Y[[i]]
    }
    gamma <- if (ns > 0) sum(S[[ns]] * Y[[ns]]) / sum(Y[[ns]] * Y[[ns]]) else 1
    r <- gamma * q
    for (i in seq_len(ns)) {
      beta <- rho[i] * sum(Y[[i]] * r); r <- r + S[[i]] * (alphas[i] - beta)
    }
    d <- -r
    # --- recherche linéaire de Wolfe (si f fourni), sinon pas unité ---
    alpha <- if (!is.null(f)) .wolfe_ls(f, grad, x, d, g, f(x)) else 1
    x_new <- x + alpha * d; g_new <- grad(x_new)
    s <- x_new - x; yv <- g_new - g
    if (sum(yv * s) > 1e-10) {                              # courbure positive
      S <- c(S, list(s)); Y <- c(Y, list(yv))
      if (length(S) > m) { S <- S[-1]; Y <- Y[-1] }
    }
    x <- x_new; g <- g_new
  }
  list(par = x, iter = k, grad_norm = sqrt(sum(g^2)),
       value = if (!is.null(f)) f(x) else NA_real_)
}
