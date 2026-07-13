# =============================================================================
# Module 18 — Méthode des moments généralisée (GMM)
# Implémente les équations de derivations/18_gmm.qmd. R base + Module 0.
# =============================================================================

#' GMM linéaire (moments d'instruments) : 2SLS et GMM efficace à deux étapes
#'
#' Moments \eqn{g_i(\beta)=Z_i(y_i-X_i^\top\beta)}. Avec `twostep = FALSE` et
#' \eqn{W=(Z^TZ)^{-1}}, renvoie le 2SLS (Prop. 18.2). Avec `twostep = TRUE`,
#' renvoie la GMM efficace (pondération \eqn{\hat S^{-1}}, robuste à
#' l'hétéroscédasticité) et le test de suridentification J (éq. 18.3).
#'
#' @param y réponse.
#' @param X régresseurs (matrice n x k, constante incluse).
#' @param Z instruments (matrice n x m, m >= k).
#' @param twostep GMM efficace à deux étapes (défaut TRUE) ; sinon 2SLS.
#' @return liste : `coefficients`, `vcov`, `se`, `J`, `J_df`, `J_pvalue`, `step`.
#' @export
gmm_linear <- function(y, X, Z, twostep = TRUE) {
  X <- as.matrix(X); Z <- as.matrix(Z); y <- as.numeric(y)
  n <- nrow(X); k <- ncol(X); m <- ncol(Z)
  if (m < k) stop("Condition d'ordre : ncol(Z) >= ncol(X).")
  A <- crossprod(Z, X)        # Z'X (m x k)
  c <- crossprod(Z, y)         # Z'y (m)
  beta_W <- function(W) as.numeric(solve(t(A) %*% W %*% A, t(A) %*% W %*% c))

  W1 <- solve(crossprod(Z))                       # (Z'Z)^{-1}
  b1 <- beta_W(W1)                                 # = 2SLS
  e1 <- as.numeric(y - X %*% b1)
  Omega <- crossprod(Z, e1^2 * Z)                  # Σ e_i^2 Z_i Z_i'

  if (twostep) {
    W2 <- solve(Omega)
    beta <- beta_W(W2)
    vcov <- solve(t(A) %*% W2 %*% A)               # (A' Omega^{-1} A)^{-1}
    r <- as.numeric(c - A %*% beta)                # moments empiriques (x n)
    J <- as.numeric(t(r) %*% W2 %*% r)             # éq. 18.3
    step <- 2L
  } else {
    beta <- b1
    # variance robuste (sandwich) du 2SLS : (A'W1A)^{-1} A'W1 Omega W1 A (A'W1A)^{-1}
    bread <- solve(t(A) %*% W1 %*% A)
    vcov <- bread %*% (t(A) %*% W1 %*% Omega %*% W1 %*% A) %*% bread
    J <- NA_real_; step <- 1L
  }
  names(beta) <- colnames(X); dimnames(vcov) <- list(colnames(X), colnames(X))
  list(coefficients = beta, vcov = vcov, se = sqrt(diag(vcov)),
       J = J, J_df = m - k, J_pvalue = if (is.na(J)) NA else pchisq(J, m - k, lower.tail = FALSE),
       step = step)
}

#' GMM générique (moments non linéaires) par optimisation
#'
#' Minimise \eqn{\bar g(\theta)^\top W\,\bar g(\theta)} (éq. 18.1). Deux étapes
#' pour la GMM efficace (\eqn{W=\hat S^{-1}}).
#'
#' @param g_fn fonction `(theta, data) -> matrice n x m` des contributions de moment.
#' @param theta0 valeur initiale.
#' @param data données passées à `g_fn`.
#' @param W matrice de pondération initiale (défaut identité).
#' @param twostep GMM efficace à deux étapes (défaut TRUE).
#' @return liste : `coefficients`, `vcov`, `se`, `J`, `J_df`, `J_pvalue`.
#' @export
gmm_fit <- function(g_fn, theta0, data, W = NULL, twostep = TRUE) {
  n <- nrow(g_fn(theta0, data)); m <- ncol(g_fn(theta0, data)); k <- length(theta0)
  gbar <- function(th) colMeans(g_fn(th, data))
  obj <- function(th, W) { gb <- gbar(th); as.numeric(t(gb) %*% W %*% gb) }
  if (is.null(W)) W <- diag(m)
  th <- optim(theta0, obj, W = W, method = "BFGS")$par
  if (twostep) {
    G1 <- g_fn(th, data); S <- crossprod(G1) / n
    W <- solve(S)
    th <- optim(th, obj, W = W, method = "BFGS")$par
  }
  # Jacobienne numérique G = d gbar / d theta (m x k)
  eps <- 1e-6
  Gmat <- vapply(seq_len(k), function(j) {
    tp <- th; tp[j] <- th[j] + eps; tm <- th; tm[j] <- th[j] - eps
    (gbar(tp) - gbar(tm)) / (2 * eps)
  }, numeric(m))
  Gmat <- matrix(Gmat, m, k)
  S <- crossprod(g_fn(th, data)) / n
  vcov <- solve(t(Gmat) %*% solve(S) %*% Gmat) / n       # (G'S^{-1}G)^{-1}/n
  gb <- gbar(th); J <- n * as.numeric(t(gb) %*% solve(S) %*% gb)
  list(coefficients = th, vcov = vcov, se = sqrt(diag(vcov)),
       J = J, J_df = m - k, J_pvalue = pchisq(J, m - k, lower.tail = FALSE))
}
