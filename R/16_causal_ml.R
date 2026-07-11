# =============================================================================
# Module 16 — Pont ML / causalité : DML et forêts causales
# Implémente les équations de derivations/16_causal_ml.qmd.
# Nuisances via les modèles des Modules 9 (forêt) et 10 (boosting), ou lm.
# =============================================================================

# Ajuste E[.y|X] sur (Xtr, ytr) et prédit sur Xte, selon la méthode de nuisance.
.fit_predict <- function(method, Xtr, ytr, Xte, B = 100L, ...) {
  Xtr <- as.data.frame(Xtr); Xte <- as.data.frame(Xte)
  df <- data.frame(.y = ytr, Xtr)
  form <- as.formula(paste(".y ~", paste(names(Xtr), collapse = " + ")))
  switch(method,
    forest = predict_forest(random_forest_fit(form, df, "anova", B = B, ...), Xte),
    boost  = predict_boost(gradient_boost(form, df, "l2", ...), Xte),
    lm = {
      Xm <- cbind(1, as.matrix(Xtr)); b <- solve_ls_qr(Xm, ytr)$coefficients
      as.numeric(cbind(1, as.matrix(Xte)) %*% b)
    },
    stop("nuisance inconnue"))
}

#' Double/Debiased ML pour le modèle partiellement linéaire (éq. 16.4)
#'
#' Résidualise Y et D par rapport à X (nuisances ML), puis estime
#' \eqn{\hat\theta = \sum \tilde D\tilde Y / \sum\tilde D^2} (score orthogonal de
#' Neyman). Cross-fitting en K blocs pour débiaiser le sur-ajustement.
#'
#' @param y réponse.
#' @param d traitement (0/1 ou continu).
#' @param X data.frame (ou matrice) de covariables.
#' @param K nombre de blocs de cross-fitting.
#' @param nuisance "forest" (M9), "boost" (M10) ou "lm".
#' @param crossfit TRUE (cross-fitting) ou FALSE (nuisances sur toutes les données).
#' @param seed graine.
#' @param ... passé aux modèles de nuisance (p.ex. B pour la forêt).
#' @return liste : `theta`, `se`, `ci`, `Ytil`, `Dtil`.
dml_plr <- function(y, d, X, K = 5L, nuisance = "forest", crossfit = TRUE,
                    seed = NULL, ...) {
  y <- as.numeric(y); d <- as.numeric(d); X <- as.data.frame(X); n <- length(y)
  if (!is.null(seed)) set.seed(seed)
  Ytil <- numeric(n); Dtil <- numeric(n)
  if (crossfit) {
    folds <- sample(rep_len(seq_len(K), n))
    for (k in seq_len(K)) {
      te <- which(folds == k); tr <- which(folds != k)
      lhat <- .fit_predict(nuisance, X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE], ...)
      mhat <- .fit_predict(nuisance, X[tr, , drop = FALSE], d[tr], X[te, , drop = FALSE], ...)
      Ytil[te] <- y[te] - lhat; Dtil[te] <- d[te] - mhat
    }
  } else {
    Ytil <- y - .fit_predict(nuisance, X, y, X, ...)
    Dtil <- d - .fit_predict(nuisance, X, d, X, ...)
  }
  theta <- sum(Dtil * Ytil) / sum(Dtil^2)                  # éq. 16.4
  psi <- Dtil * (Ytil - theta * Dtil)                       # score orthogonal
  se <- sqrt(mean(psi^2) / (n * mean(Dtil^2)^2))            # erreur standard DML
  list(theta = theta, se = se, ci = theta + c(-1, 1) * 1.96 * se,
       Ytil = Ytil, Dtil = Dtil)
}

#' T-learner : CATE par deux modèles séparés
#'
#' \eqn{\hat\tau(x) = \hat\mu_1(x) - \hat\mu_0(x)}, chaque \eqn{\mu} ajusté sur le
#' sous-groupe traité / contrôle (forêts du Module 9 par défaut).
#'
#' @param X data.frame de covariables (apprentissage).
#' @param y réponse.
#' @param d traitement (0/1).
#' @param newX covariables où prédire le CATE (défaut : X).
#' @param B nombre d'arbres.
#' @return vecteur des CATE estimés en newX.
t_learner <- function(X, y, d, newX = X, B = 200L) {
  X <- as.data.frame(X); newX <- as.data.frame(newX)
  mu1 <- .fit_predict("forest", X[d == 1, , drop = FALSE], y[d == 1], newX, B = B)
  mu0 <- .fit_predict("forest", X[d == 0, , drop = FALSE], y[d == 0], newX, B = B)
  mu1 - mu0
}

#' Arbre causal minimal (honnête)
#'
#' Splits maximisant l'hétérogénéité du traitement \eqn{n_L n_R(\hat\tau_L-
#' \hat\tau_R)^2}. Honnêteté : la structure est apprise sur une moitié, les effets
#' de feuille estimés sur l'autre. Version SIMPLIFIÉE (cf. dérivation §16.5).
#'
#' @param X data.frame de covariables.
#' @param y réponse.
#' @param d traitement (0/1).
#' @param max_depth profondeur maximale.
#' @param min_leaf effectif minimal (par groupe de traitement) par feuille.
#' @param seed graine (partage honnête).
#' @return objet `causal_tree`.
causal_tree <- function(X, y, d, max_depth = 3L, min_leaf = 10L, seed = NULL) {
  X <- as.matrix(X); y <- as.numeric(y); d <- as.numeric(d); n <- nrow(X)
  if (!is.null(seed)) set.seed(seed)
  struct <- sample.int(n, n %/% 2)                          # échantillon de structure
  est <- setdiff(seq_len(n), struct)                        # échantillon d'estimation
  vars <- colnames(X); if (is.null(vars)) vars <- paste0("V", seq_len(ncol(X)))

  tau_of <- function(idx) {
    yt <- y[idx][d[idx] == 1]; yc <- y[idx][d[idx] == 0]
    if (length(yt) == 0 || length(yc) == 0) return(NA_real_)
    mean(yt) - mean(yc)
  }
  enough <- function(idx) sum(d[idx] == 1) >= min_leaf && sum(d[idx] == 0) >= min_leaf

  build <- function(idx_s, depth) {
    leaf <- list(leaf = TRUE)
    if (depth >= max_depth || !enough(idx_s)) return(leaf)
    best <- list(gain = 0)
    for (j in seq_len(ncol(X))) {
      xs <- sort(unique(X[idx_s, j]))
      if (length(xs) < 2) next
      for (thr in (xs[-length(xs)] + xs[-1]) / 2) {
        L <- idx_s[X[idx_s, j] <= thr]; R <- idx_s[X[idx_s, j] > thr]
        if (!enough(L) || !enough(R)) next
        g <- length(L) * length(R) * (tau_of(L) - tau_of(R))^2
        if (is.finite(g) && g > best$gain) best <- list(gain = g, var = j, val = thr)
      }
    }
    if (best$gain <= 0) return(leaf)
    goL <- X[idx_s, best$var] <= best$val
    list(leaf = FALSE, var_idx = best$var, val = best$val,
         left = build(idx_s[goL], depth + 1L), right = build(idx_s[!goL], depth + 1L))
  }
  tree <- build(struct, 0L)

  # Estimation honnête des effets de feuille sur l'échantillon d'estimation.
  assign_leaf <- function(node, idx_e) {
    if (node$leaf) { node$tau <- tau_of(idx_e); return(node) }
    goL <- X[idx_e, node$var_idx] <= node$val
    node$left <- assign_leaf(node$left, idx_e[goL])
    node$right <- assign_leaf(node$right, idx_e[!goL])
    node
  }
  tree <- assign_leaf(tree, est)
  structure(list(tree = tree, vars = vars, tau_root = tau_of(est)), class = "causal_tree")
}

#' Prédiction du CATE par un arbre causal
#' @param object objet `causal_tree`. @param newdata covariables.
#' @return vecteur des CATE prédits.
predict_causal_tree <- function(object, newdata) {
  X <- as.matrix(newdata)
  descend <- function(node, xrow) {
    if (node$leaf) return(if (is.na(node$tau)) object$tau_root else node$tau)
    if (xrow[node$var_idx] <= node$val) descend(node$left, xrow) else descend(node$right, xrow)
  }
  vapply(seq_len(nrow(X)), function(i) descend(object$tree, X[i, ]), numeric(1))
}
