# =============================================================================
# Module 35 — Entrainement des reseaux : optimiseurs modernes et regularisation
# Implemente les equations de derivations/35_nn_training.qmd. R base.
# Le Module 12 entrainait un MLP par SGD nu. En pratique on utilise des
# optimiseurs ADAPTATIFS (Adam) et des couches de REGULARISATION (dropout, batch
# normalization) qui accelerent et stabilisent l'apprentissage profond.
# =============================================================================

#' Optimiseur Adam (moments adaptatifs)
#'
#' Combine momentum (moment 1) et mise a l'echelle par la variance des gradients
#' (moment 2), avec correction de biais : \eqn{m_t=\beta_1 m_{t-1}+(1-\beta_1)g_t},
#' \eqn{v_t=\beta_2 v_{t-1}+(1-\beta_2)g_t^2}, pas
#' \eqn{x_{t+1}=x_t-\eta\,\hat m_t/(\sqrt{\hat v_t}+\epsilon)}.
#'
#' @param grad fonction `x -> gradient`.
#' @param x0 point initial.
#' @param lr pas d'apprentissage.
#' @param beta1,beta2 taux de decroissance des moments.
#' @param eps stabilisateur numerique.
#' @param max_iter,tol arret.
#' @return liste : `par`, `iters`, `path` (valeurs successives).
#' @export
optim_adam <- function(grad, x0, lr = 0.01, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
                       max_iter = 5000L, tol = 1e-8) {
  x <- x0; m <- v <- rep(0, length(x0)); path <- matrix(x0, 1)
  for (t in seq_len(max_iter)) {
    g <- grad(x); m <- beta1 * m + (1 - beta1) * g; v <- beta2 * v + (1 - beta2) * g^2
    mh <- m / (1 - beta1^t); vh <- v / (1 - beta2^t)
    step <- lr * mh / (sqrt(vh) + eps); x <- x - step
    path <- rbind(path, x)
    if (max(abs(step)) < tol) break
  }
  list(par = x, iters = t, path = path)
}

#' Optimiseur RMSprop
#'
#' \eqn{v_t=\rho v_{t-1}+(1-\rho)g_t^2}, \eqn{x_{t+1}=x_t-\eta g_t/(\sqrt{v_t}+\epsilon)}.
#'
#' @param grad,x0,lr cf. `optim_adam` ; @param rho decroissance ; @param eps,max_iter,tol arret.
#' @return liste : `par`, `iters`.
#' @export
optim_rmsprop <- function(grad, x0, lr = 0.01, rho = 0.9, eps = 1e-8, max_iter = 5000L, tol = 1e-8) {
  x <- x0; v <- rep(0, length(x0))
  for (t in seq_len(max_iter)) {
    g <- grad(x); v <- rho * v + (1 - rho) * g^2
    step <- lr * g / (sqrt(v) + eps); x <- x - step
    if (max(abs(step)) < tol) break
  }
  list(par = x, iters = t)
}

#' Descente de gradient a momentum (boule pesante de Polyak)
#'
#' \eqn{u_t=\mu u_{t-1}-\eta g_t}, \eqn{x_{t+1}=x_t+u_t}.
#'
#' @param grad,x0 cf. `optim_adam` ; @param step pas ; @param momentum \eqn{\mu} ;
#' @param max_iter,tol arret.
#' @return liste : `par`, `iters`.
#' @export
optim_momentum <- function(grad, x0, step = 0.01, momentum = 0.9, max_iter = 5000L, tol = 1e-8) {
  x <- x0; u <- rep(0, length(x0))
  for (t in seq_len(max_iter)) {
    g <- grad(x); u <- momentum * u - step * g; x <- x + u
    if (max(abs(u)) < tol) break
  }
  list(par = x, iters = t)
}

#' Dropout (inverted dropout)
#'
#' A l'apprentissage, met a zero une fraction `rate` des unites et remet a
#' l'echelle par \eqn{1/(1-\text{rate})} (l'esperance est preservee). A
#' l'evaluation (`training=FALSE`), l'identite. Regularise en empechant la
#' co-adaptation des neurones.
#'
#' @param x activations (vecteur ou matrice).
#' @param rate probabilite d'extinction.
#' @param training TRUE (apprentissage) ou FALSE (evaluation).
#' @return liste : `out`, `mask`.
#' @export
dropout <- function(x, rate = 0.5, training = TRUE) {
  if (!training || rate == 0) return(list(out = x, mask = array(1, dim = dim(x) %||% length(x))))
  mask <- (matrix(runif(length(x)), nrow = NROW(x)) > rate) / (1 - rate)
  if (is.null(dim(x))) mask <- as.numeric(mask)
  list(out = x * mask, mask = mask)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Batch normalization (passe avant)
#'
#' Normalise chaque colonne (feature) sur le batch : \eqn{\hat x=(x-\mu)/\sqrt{
#' \sigma^2+\epsilon}}, puis \eqn{y=\gamma\hat x+\beta}. Stabilise et accelere
#' l'apprentissage (reduit le decalage de covariance interne).
#'
#' @param X matrice n x d (batch).
#' @param gamma,beta parametres d'echelle et de decalage (longueur d).
#' @param eps stabilisateur.
#' @return liste : `out`, `cache` (pour la retropropagation).
#' @export
batch_norm <- function(X, gamma = rep(1, ncol(X)), beta = rep(0, ncol(X)), eps = 1e-5) {
  mu <- colMeans(X); xc <- sweep(X, 2, mu); vr <- colMeans(xc^2)
  istd <- 1 / sqrt(vr + eps); xhat <- sweep(xc, 2, istd, "*")
  out <- sweep(sweep(xhat, 2, gamma, "*"), 2, beta, "+")
  list(out = out, cache = list(xhat = xhat, gamma = gamma, istd = istd, xc = xc, n = nrow(X)))
}

#' Batch normalization (passe arriere)
#'
#' Gradients de la perte par rapport a \eqn{X}, \eqn{\gamma}, \eqn{\beta} etant
#' donne `dout` (gradient en sortie).
#'
#' @param dout gradient en sortie (n x d).
#' @param cache sortie de `batch_norm`.
#' @return liste : `dX`, `dgamma`, `dbeta`.
#' @export
batch_norm_backward <- function(dout, cache) {
  n <- cache$n; xhat <- cache$xhat; istd <- cache$istd; xc <- cache$xc; gamma <- cache$gamma
  dbeta <- colSums(dout); dgamma <- colSums(dout * xhat)
  dxhat <- sweep(dout, 2, gamma, "*")
  dvar <- colSums(dxhat * xc) * -0.5 * istd^3
  dmu <- colSums(sweep(dxhat, 2, -istd, "*")) + dvar * colMeans(-2 * xc)
  dX <- sweep(dxhat, 2, istd, "*") + sweep(xc, 2, 2 * dvar / n, "*") + rep(dmu / n, each = n)
  list(dX = dX, dgamma = dgamma, dbeta = dbeta)
}
