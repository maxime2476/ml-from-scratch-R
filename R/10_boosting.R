# =============================================================================
# Module 10 — Gradient boosting
# Implémente les équations de derivations/10_boosting.qmd.
# Arbres de base = cart_fit (Module 8, method "anova") sur les pseudo-résidus ;
# raffinement des feuilles par le poids de Newton (éq. 10.7) en option.
# =============================================================================

# Raffinement des poids de feuille par le pas de Newton (éq. 10.7).
# Route chaque échantillon vers sa feuille et y pose w = -sum(g)/(sum(h)+lambda).
.refit_leaves_newton <- function(node, idx, X, g, h, lambda) {
  if (node$leaf) {
    node$pred <- if (length(idx) == 0) 0 else -sum(g[idx]) / (sum(h[idx]) + lambda)
    return(node)
  }
  go_left <- X[idx, node$var_idx] <= node$val
  node$left  <- .refit_leaves_newton(node$left,  idx[go_left],  X, g, h, lambda)
  node$right <- .refit_leaves_newton(node$right, idx[!go_left], X, g, h, lambda)
  node
}

#' Gradient boosting (descente de gradient fonctionnelle, éq. 10.1-10.2)
#'
#' Ajuste additivement M arbres aux pseudo-résidus. `loss = "l2"` (résidu, éq.
#' 10.3) ou `"logloss"` (y - p, éq. 10.5). Avec `newton = TRUE`, les poids de
#' feuille sont ceux de Newton (éq. 10.7) — pour la log-loss, la mise à jour de
#' Friedman (comme `gbm`).
#'
#' @param formula formule (prédicteurs numériques).
#' @param data data.frame.
#' @param loss "l2" (régression) ou "logloss" (classification binaire y in {0,1}).
#' @param M nombre d'arbres.
#' @param nu taux d'apprentissage (shrinkage) dans (0,1].
#' @param max_depth,min_split,min_leaf hyperparamètres des arbres de base.
#' @param lambda régularisation L2 des poids de feuille (éq. 10.7).
#' @param newton raffiner les feuilles par le pas de Newton (défaut TRUE).
#' @return objet `boost`.
gradient_boost <- function(formula, data, loss = c("l2", "logloss"), M = 100L,
                           nu = 0.1, max_depth = 3L, min_split = 10L, min_leaf = 5L,
                           lambda = 0, newton = TRUE) {
  loss <- match.arg(loss)
  mf <- model.frame(formula, data)
  y <- as.numeric(model.response(mf))
  vars <- colnames(mf)[-1]
  X <- as.matrix(mf[, -1, drop = FALSE]); storage.mode(X) <- "double"
  n <- nrow(X)
  dtr <- as.data.frame(X); names(dtr) <- vars
  tree_formula <- as.formula(paste(".resid ~", paste(vars, collapse = " + ")))

  F0 <- if (loss == "l2") mean(y) else { pb <- mean(y); log(pb / (1 - pb)) }
  Fcur <- rep(F0, n)
  trees <- vector("list", M)

  for (m in seq_len(M)) {
    if (loss == "l2") { g <- -(y - Fcur); h <- rep(1, n) }
    else { p <- 1 / (1 + exp(-Fcur)); g <- p - y; h <- p * (1 - p) }
    dtr$.resid <- -g                                # pseudo-résidu (éq. 10.1)
    tr <- cart_fit(tree_formula, dtr, "anova", max_depth = max_depth,
                   min_split = min_split, min_leaf = min_leaf)
    if (newton) tr$tree <- .refit_leaves_newton(tr$tree, seq_len(n), X, g, h, lambda)
    Fcur <- Fcur + nu * predict_cart(tr, dtr)       # mise à jour (éq. 10.2)
    trees[[m]] <- tr
  }
  structure(list(trees = trees, F0 = F0, nu = nu, loss = loss, vars = vars, M = M),
            class = "boost")
}

#' Prédiction d'un modèle de boosting
#'
#' @param object objet `boost`.
#' @param newdata data.frame des prédicteurs.
#' @param type "response" (valeur / probabilité), "link" (score F) ou "class"
#'   (0/1, log-loss uniquement).
#' @param n_trees nombre d'arbres à utiliser (défaut : tous).
#' @return vecteur de prédictions.
predict_boost <- function(object, newdata, type = c("response", "link", "class"),
                          n_trees = object$M) {
  type <- match.arg(type)
  F <- rep(object$F0, nrow(newdata))
  for (m in seq_len(n_trees)) F <- F + object$nu * predict_cart(object$trees[[m]], newdata)
  if (type == "link") return(F)
  if (object$loss == "l2") return(F)
  p <- 1 / (1 + exp(-F))
  if (type == "class") as.integer(p >= 0.5) else p
}

#' Trajectoire de la perte d'entraînement/test selon le nombre d'arbres
#'
#' @param object objet `boost`.
#' @param data data.frame des prédicteurs.
#' @param y réponse correspondante.
#' @return vecteur de perte (moyenne de \eqn{\ell}) après m arbres, m = 1..M.
boost_loss_path <- function(object, data, y) {
  y <- as.numeric(y); n <- nrow(data)
  F <- rep(object$F0, n); losses <- numeric(object$M)
  for (m in seq_len(object$M)) {
    F <- F + object$nu * predict_cart(object$trees[[m]], data)
    losses[m] <- if (object$loss == "l2") mean((y - F)^2)
                 else mean(log(1 + exp(F)) - y * F)   # log-loss moyenne
  }
  losses
}
