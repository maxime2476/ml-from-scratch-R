# =============================================================================
# Module 12 — Perceptron multicouche (une couche cachée)
# Implémente les équations de derivations/12_mlp.qmd. R base, from scratch.
# =============================================================================

# ---- Activations et dérivées ------------------------------------------------
.act <- function(Z, kind) switch(kind,
  tanh = tanh(Z), relu = pmax(Z, 0), sigmoid = 1 / (1 + exp(-Z)),
  stop("activation inconnue"))
.act_grad <- function(Z, kind) switch(kind,
  tanh = 1 - tanh(Z)^2, relu = (Z > 0) * 1,
  sigmoid = { s <- 1 / (1 + exp(-Z)); s * (1 - s) })
.output <- function(Z2, loss) if (loss == "mse") Z2 else 1 / (1 + exp(-Z2))

#' Passe avant du MLP (éq. 12.1)
#'
#' @param params liste `W1, b1, W2, b2`.
#' @param X matrice n x d0.
#' @param activation activation de la couche cachée ("tanh", "relu", "sigmoid").
#' @return liste `Z1, A1, Z2`.
mlp_forward <- function(params, X, activation) {
  Z1 <- sweep(X %*% params$W1, 2, params$b1, "+")
  A1 <- .act(Z1, activation)
  Z2 <- sweep(A1 %*% params$W2, 2, params$b2, "+")
  list(Z1 = Z1, A1 = A1, Z2 = Z2)
}

#' Perte moyenne du MLP
#'
#' @param Yhat prédictions (sortie du réseau).
#' @param Y cible.
#' @param loss "mse" (régression) ou "logloss" (classification binaire).
#' @return la perte moyenne.
mlp_loss <- function(Yhat, Y, loss) {
  n <- nrow(Y)
  if (loss == "mse") sum((Yhat - Y)^2) / (2 * n)
  else { eps <- 1e-12; -sum(Y * log(Yhat + eps) + (1 - Y) * log(1 - Yhat + eps)) / n }
}

#' Rétropropagation : gradients analytiques (éq. 12.2-12.5)
#'
#' @param params liste des paramètres.
#' @param X matrice n x d0.
#' @param Y cible n x d2.
#' @param activation activation cachée.
#' @param loss "mse" ou "logloss".
#' @return liste des gradients `W1, b1, W2, b2`.
mlp_backward <- function(params, X, Y, activation, loss) {
  n <- nrow(X)
  cache <- mlp_forward(params, X, activation)
  Yhat <- .output(cache$Z2, loss)
  d2 <- (Yhat - Y) / n                                   # éq. 12.2
  gW2 <- crossprod(cache$A1, d2)                          # A1' d2 (éq. 12.3)
  gb2 <- colSums(d2)
  d1 <- (d2 %*% t(params$W2)) * .act_grad(cache$Z1, activation)   # éq. 12.4
  gW1 <- crossprod(X, d1)                                 # X' d1 (éq. 12.5)
  gb1 <- colSums(d1)
  list(W1 = gW1, b1 = gb1, W2 = gW2, b2 = gb2)
}

#' Entraînement d'un MLP par SGD (éq. 12.6)
#'
#' @param X matrice n x d0.
#' @param y cible (vecteur ou matrice n x d2).
#' @param hidden largeur de la couche cachée d1.
#' @param activation "tanh" (défaut), "relu" ou "sigmoid".
#' @param loss "mse" ou "logloss".
#' @param epochs nombre d'époques.
#' @param lr taux d'apprentissage.
#' @param batch taille de mini-lot.
#' @param seed graine (initialisation + mélange).
#' @return objet `mlp` : `params`, `activation`, `loss`, `loss_hist`, `d`.
mlp_fit <- function(X, y, hidden = 8L, activation = "tanh",
                    loss = c("mse", "logloss"), epochs = 200L, lr = 0.05,
                    batch = 32L, seed = NULL) {
  loss <- match.arg(loss)
  X <- as.matrix(X); Y <- as.matrix(y)
  n <- nrow(X); d0 <- ncol(X); d1 <- hidden; d2 <- ncol(Y)
  if (!is.null(seed)) set.seed(seed)
  params <- list(
    W1 = matrix(rnorm(d0 * d1) * sqrt(1 / d0), d0, d1), b1 = rep(0, d1),
    W2 = matrix(rnorm(d1 * d2) * sqrt(1 / d1), d1, d2), b2 = rep(0, d2))
  loss_hist <- numeric(epochs)
  for (e in seq_len(epochs)) {
    idx <- sample.int(n)
    for (s in seq(1L, n, by = batch)) {
      bi <- idx[s:min(s + batch - 1L, n)]
      g <- mlp_backward(params, X[bi, , drop = FALSE], Y[bi, , drop = FALSE], activation, loss)
      params$W1 <- params$W1 - lr * g$W1; params$b1 <- params$b1 - lr * g$b1
      params$W2 <- params$W2 - lr * g$W2; params$b2 <- params$b2 - lr * g$b2
    }
    Yhat <- .output(mlp_forward(params, X, activation)$Z2, loss)
    loss_hist[e] <- mlp_loss(Yhat, Y, loss)
  }
  structure(list(params = params, activation = activation, loss = loss,
                 loss_hist = loss_hist, d = c(d0, d1, d2)), class = "mlp")
}

#' Prédiction d'un MLP
#'
#' @param model objet `mlp`.
#' @param newdata matrice des prédicteurs.
#' @param type "response" (valeur/probabilité) ou "class" (0/1, log-loss).
#' @return vecteur/matrice de prédictions.
predict_mlp <- function(model, newdata, type = c("response", "class")) {
  type <- match.arg(type)
  Z2 <- mlp_forward(model$params, as.matrix(newdata), model$activation)$Z2
  out <- .output(Z2, model$loss)
  if (model$loss == "logloss" && type == "class") (out >= 0.5) * 1 else drop(out)
}

#' Gradient numérique par différences finies centrées (éq. 12.7)
#'
#' Sert à VÉRIFIER la rétropropagation. Perturbe chaque paramètre de +/- eps et
#' approche le gradient par \eqn{(L(\theta+\varepsilon)-L(\theta-\varepsilon))/(2\varepsilon)}.
#'
#' @param params liste des paramètres.
#' @param X,Y données.
#' @param activation activation cachée.
#' @param loss "mse" ou "logloss".
#' @param eps pas de différence finie.
#' @return liste des gradients numériques (même structure que `params`).
mlp_numgrad <- function(params, X, Y, activation, loss, eps = 1e-6) {
  Lf <- function(p) mlp_loss(.output(mlp_forward(p, X, activation)$Z2, loss), Y, loss)
  out <- params
  for (nm in names(params)) {
    v <- params[[nm]]; g <- numeric(length(v))
    for (k in seq_along(v)) {
      pp <- params; pp[[nm]][k] <- v[k] + eps; Lp <- Lf(pp)
      pp[[nm]][k] <- v[k] - eps; Lm <- Lf(pp)
      g[k] <- (Lp - Lm) / (2 * eps)
    }
    out[[nm]] <- if (is.matrix(v)) matrix(g, nrow(v), ncol(v)) else g
  }
  out
}
