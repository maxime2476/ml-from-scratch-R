# =============================================================================
# Module 8 — Arbres de décision (CART)
# Implémente les équations de derivations/08_cart.qmd. R base uniquement.
# Prédicteurs numériques ; classification (Gini/entropie) et régression (variance).
# =============================================================================

# ---- Mesures d'impureté (éq. 8.1-8.2) ---------------------------------------

#' Impureté de Gini d'un vecteur d'étiquettes (éq. 8.1)
#' @param y étiquettes (facteur/vecteur). @return \eqn{1-\sum_c p_c^2}.
impurity_gini <- function(y) { p <- prop.table(table(y)); 1 - sum(p^2) }

#' Entropie d'un vecteur d'étiquettes (éq. 8.1)
#' @param y étiquettes. @return \eqn{-\sum_c p_c\log p_c}.
impurity_entropy <- function(y) { p <- prop.table(table(y)); p <- p[p > 0]; -sum(p * log(p)) }

#' Impureté de variance (régression, éq. 8.2)
#' @param y réponse numérique. @return variance de population intra-nœud.
impurity_variance <- function(y) mean((y - mean(y))^2)

.imp_counts <- function(cnt, kind) {
  n <- sum(cnt); if (n == 0) return(0)
  p <- cnt / n; p <- p[p > 0]
  if (kind == "gini") 1 - sum(p^2) else -sum(p * log(p))
}

# ---- Recherche du meilleur split (éq. 8.3) ----------------------------------

#' Meilleur split d'un nœud (éq. 8.3)
#'
#' Balayage incrémental sur chaque variable ; renvoie le couple (variable, seuil)
#' minimisant l'impureté pondérée des enfants (= maximisant \eqn{\Delta I}).
#'
#' @param X matrice de prédicteurs (numériques) du nœud.
#' @param y réponse du nœud.
#' @param method "class" ou "anova".
#' @param kind pour la classification : "gini" ou "entropy".
#' @param min_leaf effectif minimal par feuille.
#' @param classes niveaux de classe (classification).
#' @return liste `gain`, `var` (indice), `val` (seuil) ; `gain = -Inf` si aucun split.
best_split <- function(X, y, method, kind = "gini", min_leaf = 1L, classes = NULL) {
  n <- nrow(X); best <- list(gain = -Inf, var = NA_integer_, val = NA_real_)
  if (method == "class") {
    yc <- match(y, classes); C <- length(classes)
    tot <- tabulate(yc, C); parent <- .imp_counts(tot, kind)
    for (j in seq_len(ncol(X))) {
      ord <- order(X[, j]); xs <- X[ord, j]; yy <- yc[ord]
      left <- integer(C)
      for (i in seq_len(n - 1L)) {
        left[yy[i]] <- left[yy[i]] + 1L
        if (xs[i] == xs[i + 1L]) next
        nl <- i; nr <- n - i
        if (nl < min_leaf || nr < min_leaf) next
        imp <- (nl * .imp_counts(left, kind) + nr * .imp_counts(tot - left, kind)) / n
        gain <- parent - imp
        if (gain > best$gain) best <- list(gain = gain, var = j, val = unname((xs[i] + xs[i + 1L]) / 2))
      }
    }
  } else {
    tsum <- sum(y); tsq <- sum(y^2); sse_tot <- tsq - tsum^2 / n
    for (j in seq_len(ncol(X))) {
      ord <- order(X[, j]); xs <- X[ord, j]; ys <- y[ord]
      sl <- 0; sql <- 0
      for (i in seq_len(n - 1L)) {
        sl <- sl + ys[i]; sql <- sql + ys[i]^2
        if (xs[i] == xs[i + 1L]) next
        nl <- i; nr <- n - i
        if (nl < min_leaf || nr < min_leaf) next
        sse_l <- sql - sl^2 / nl
        sse_r <- (tsq - sql) - (tsum - sl)^2 / nr
        gain <- (sse_tot - sse_l - sse_r) / n
        if (gain > best$gain) best <- list(gain = gain, var = j, val = unname((xs[i] + xs[i + 1L]) / 2))
      }
    }
  }
  best
}

# ---- Croissance récursive ----------------------------------------------------

#' Ajustement d'un arbre CART
#'
#' Croissance récursive gloutonne (§8.3-8.4) : à chaque nœud, `best_split`
#' choisit le split optimal ; arrêt sur profondeur, effectif, pureté ou gain nul.
#'
#' @param formula formule façon `rpart` (prédicteurs numériques).
#' @param data data.frame.
#' @param method "class" (Gini/entropie) ou "anova" (variance).
#' @param kind impureté de classification : "gini" (défaut) ou "entropy".
#' @param max_depth profondeur maximale.
#' @param min_split effectif minimal pour tenter un split.
#' @param min_leaf effectif minimal par feuille.
#' @param min_gain gain d'impureté minimal pour accepter un split.
#' @return objet `cart` (arbre + métadonnées).
cart_fit <- function(formula, data, method = c("class", "anova"), kind = "gini",
                     max_depth = 30L, min_split = 20L, min_leaf = 7L, min_gain = 1e-9) {
  method <- match.arg(method)
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- as.matrix(mf[, -1, drop = FALSE])
  storage.mode(X) <- "double"
  vars <- colnames(X)
  classes <- if (method == "class") levels(as.factor(y)) else NULL
  N <- nrow(X)

  node_pred <- function(yn) {
    if (method == "class") {
      tb <- table(factor(yn, levels = classes))
      list(pred = names(tb)[which.max(tb)], prob = as.numeric(prop.table(tb)))
    } else list(pred = mean(yn), prob = NULL)
  }
  node_imp <- function(yn) {
    if (method == "class") (if (kind == "gini") impurity_gini(yn) else impurity_entropy(yn))
    else impurity_variance(yn)
  }

  build <- function(idx, depth) {
    yn <- y[idx]; np <- node_pred(yn); imp <- node_imp(yn)
    leaf <- list(leaf = TRUE, pred = np$pred, prob = np$prob, n = length(idx), impurity = imp)
    if (depth >= max_depth || length(idx) < min_split || imp <= 0) return(leaf)
    bs <- best_split(X[idx, , drop = FALSE], yn, method, kind, min_leaf, classes)
    if (!is.finite(bs$gain) || bs$gain <= min_gain) return(leaf)
    go_left <- X[idx, bs$var] <= bs$val
    list(leaf = FALSE, var = vars[bs$var], var_idx = bs$var, val = bs$val,
         gain = bs$gain, n = length(idx), impurity = imp,
         pred = np$pred, prob = np$prob,        # prédiction si on élague ici
         left = build(idx[go_left], depth + 1L), right = build(idx[!go_left], depth + 1L))
  }

  structure(list(tree = build(seq_len(N), 0L), method = method, kind = kind,
                 vars = vars, classes = classes, N = N), class = "cart")
}

#' Prédiction d'un arbre CART
#'
#' @param object objet `cart`.
#' @param newdata data.frame contenant les prédicteurs.
#' @return vecteur de prédictions (classes pour "class", moyennes pour "anova").
predict_cart <- function(object, newdata) {
  X <- as.matrix(newdata[, object$vars, drop = FALSE]); storage.mode(X) <- "double"
  descend <- function(node, xrow) {
    if (node$leaf) return(node$pred)
    if (xrow[node$var_idx] <= node$val) descend(node$left, xrow) else descend(node$right, xrow)
  }
  out <- vapply(seq_len(nrow(X)), function(i) as.character(descend(object$tree, X[i, ])), character(1))
  if (object$method == "class") factor(out, levels = object$classes) else as.numeric(out)
}

#' Élagage coût-complexité (éq. 8.4-8.5)
#'
#' Effondre récursivement les sous-arbres dont le maintien n'améliore pas
#' \eqn{R_\alpha(T)=R(T)+\alpha|T|} (weakest-link, Prop. 8.2).
#'
#' @param object objet `cart`.
#' @param alpha coût par feuille \eqn{\alpha \ge 0}.
#' @return objet `cart` élagué.
cost_complexity_prune <- function(object, alpha) {
  N <- object$N
  prune <- function(node) {
    if (node$leaf) return(list(node = node, R = (node$n / N) * node$impurity, leaves = 1L))
    L <- prune(node$left); Rr <- prune(node$right)
    R_sub <- L$R + Rr$R; leaves <- L$leaves + Rr$leaves
    R_collapse <- (node$n / N) * node$impurity + alpha              # feuille + alpha*1
    R_keep <- R_sub + alpha * leaves
    if (R_collapse <= R_keep) {                                     # effondrement favorable
      list(node = list(leaf = TRUE, pred = node$pred, prob = node$prob,
                       n = node$n, impurity = node$impurity),
           R = (node$n / N) * node$impurity, leaves = 1L)
    } else {
      node$left <- L$node; node$right <- Rr$node
      list(node = node, R = R_sub, leaves = leaves)
    }
  }
  object$tree <- prune(object$tree)$node
  object
}

#' Nombre de feuilles d'un arbre CART
#' @param object objet `cart`. @return nombre de feuilles.
n_leaves <- function(object) {
  count <- function(node) if (node$leaf) 1L else count(node$left) + count(node$right)
  count(object$tree)
}
