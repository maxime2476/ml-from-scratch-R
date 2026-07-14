# =============================================================================
# Module 41 — Clustering avance : hierarchique, DBSCAN, spectral
# Implemente les equations de derivations/41_clustering.qmd. R base.
# Le k-means (Module 11) suppose des groupes convexes et spheriques. Trois
# approches le depassent : la HIERARCHIE (arbre de fusions), la DENSITE (DBSCAN,
# formes arbitraires + bruit) et le SPECTRE du graphe (clusters non convexes).
# =============================================================================

#' Clustering hierarchique agglomeratif (assignation a k groupes)
#'
#' Fusionne iterativement les deux groupes les plus proches selon la **liaison**
#' (mise a jour de Lance-Williams), jusqu'a `k` groupes. Liaisons : "complete"
#' (diametre), "single" (chaine), "average" (moyenne).
#'
#' @param X matrice n x p (ou objet `dist`) ; @param k nombre de groupes ;
#' @param linkage "complete", "single" ou "average".
#' @return vecteur d'assignation (longueur n).
#' @export
agglomerative <- function(X, k, linkage = c("complete", "single", "average")) {
  linkage <- match.arg(linkage)
  D <- if (inherits(X, "dist")) X else dist(X)
  n <- attr(D, "Size"); Dm <- as.matrix(D); diag(Dm) <- Inf
  members <- as.list(seq_len(n)); sizes <- rep(1, n); active <- rep(TRUE, n)
  while (sum(active) > k) {
    act <- which(active); sub <- Dm[act, act]
    wm <- which(sub == min(sub), arr.ind = TRUE)[1, ]; i <- act[wm[1]]; j <- act[wm[2]]
    for (m in act) if (m != i && m != j) {
      Dm[i, m] <- Dm[m, i] <- switch(linkage,
        complete = max(Dm[i, m], Dm[j, m]), single = min(Dm[i, m], Dm[j, m]),
        average = (sizes[i] * Dm[i, m] + sizes[j] * Dm[j, m]) / (sizes[i] + sizes[j]))
    }
    members[[i]] <- c(members[[i]], members[[j]]); sizes[i] <- sizes[i] + sizes[j]; active[j] <- FALSE
  }
  assign <- integer(n); rt <- which(active)
  for (c in seq_along(rt)) assign[members[[rt[c]]]] <- c
  assign
}

#' DBSCAN (clustering base sur la densite)
#'
#' Un point est **coeur** s'il a au moins `minPts` voisins dans un rayon `eps`.
#' Les clusters sont des composantes connexes de points-coeur (et leurs voisins) ;
#' les points isoles sont du **bruit** (label 0). Detecte des formes arbitraires,
#' sans fixer le nombre de clusters.
#'
#' @param X matrice n x p ; @param eps rayon ; @param minPts densite minimale.
#' @return vecteur d'etiquettes (0 = bruit).
#' @export
dbscan_fit <- function(X, eps, minPts = 5L) {
  X <- as.matrix(X); n <- nrow(X); Dm <- as.matrix(dist(X))
  lab <- rep(0L, n); visited <- logical(n); cid <- 0L
  neigh <- function(i) which(Dm[i, ] <= eps)
  for (i in seq_len(n)) {
    if (visited[i]) next; visited[i] <- TRUE; N <- neigh(i)
    if (length(N) < minPts) next                            # reste bruit (0) pour l'instant
    cid <- cid + 1L; lab[i] <- cid; seeds <- setdiff(N, i)
    while (length(seeds) > 0) {
      j <- seeds[1]; seeds <- seeds[-1]
      if (!visited[j]) { visited[j] <- TRUE; Nj <- neigh(j)
        if (length(Nj) >= minPts) seeds <- union(seeds, setdiff(Nj, which(visited))) }
      if (lab[j] == 0L) lab[j] <- cid
    }
  }
  lab
}

#' Clustering spectral
#'
#' Construit un graphe de similarite (noyau RBF), calcule le **Laplacien
#' normalise** \eqn{L=I-D^{-1/2}WD^{-1/2}}, plonge les points dans les `k` premiers
#' vecteurs propres, puis applique k-means. Separe des clusters **non convexes**
#' (spirales, cercles) que le k-means echoue a distinguer.
#'
#' @param X matrice n x p ; @param k nombre de clusters ; @param gamma echelle RBF.
#' @return vecteur d'assignation.
#' @export
spectral_clustering <- function(X, k, gamma = 1) {
  X <- as.matrix(X); n <- nrow(X)
  D2 <- outer(rowSums(X^2), rowSums(X^2), "+") - 2 * X %*% t(X)
  W <- exp(-gamma * pmax(D2, 0)); diag(W) <- 0
  d <- rowSums(W); Dm12 <- 1 / sqrt(d)
  L <- diag(n) - (Dm12 * W) * rep(Dm12, each = n)           # Laplacien normalise symetrique
  ev <- eigen(L, symmetric = TRUE)
  U <- ev$vectors[, (n - k + 1):n, drop = FALSE]            # k plus petits (vp les plus faibles)
  U <- U / sqrt(rowSums(U^2) + 1e-12)                       # normalisation par ligne (Ng-Jordan-Weiss)
  km <- kmeans(U, centers = k, nstart = 10)
  km$cluster
}
