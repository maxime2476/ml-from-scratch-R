# =============================================================================
# Module 22 — Lasso débiaisé (désparsifié)
# Implémente les équations de derivations/22_debiased_lasso.qmd.
# Réutilise lasso_fit (Module 4) et la standardisation. R base + Module 4.
# =============================================================================

# lambda théorique (convention lasso_fit : (1/2)||.||^2 + lambda||.||_1 sur X
# standardisée, ||x_j||^2 = n) : lambda ~ c * sigma * sqrt(2 n log p).
.lambda_theo <- function(n, p, sigma, c = 1.1) c * sigma * sqrt(2 * n * log(p))

#' Lasso débiaisé / désparsifié (inférence haute dimension valide)
#'
#' Corrige le biais de rétrécissement du lasso par une **projection de faible
#' dimension** (Zhang-Zhang 2014) : pour chaque coordonnée j, on partiale les
#' autres covariables hors de \eqn{x_j} par un lasso « nodewise », et on utilise
#' le résidu \eqn{\hat\tau_j} comme score orthogonal (même idée que Neyman, M16) :
#' \eqn{\hat\beta^d_j=\hat\beta_j + \hat\tau_j^\top(y-X\hat\beta)/(\hat\tau_j^\top x_j)}.
#' Fournit des intervalles de confiance **valides** même si p > n — là où le
#' t-test post-lasso naïf échoue (Module 14).
#'
#' @param X matrice n x p.
#' @param y réponse.
#' @param lambda pénalité du lasso principal (défaut : théorique).
#' @param lambda_node pénalité des lasso nodewise (défaut : théorique).
#' @param targets indices des coordonnées à débiaiser (défaut : toutes).
#' @param sigma écart-type du bruit (défaut : estimé sur les résidus du lasso).
#' @return liste : `estimate`, `se`, `lower`, `upper`, `beta_lasso`, `sigma`.
#' @export
debiased_lasso <- function(X, y, lambda = NULL, lambda_node = NULL,
                           targets = NULL, sigma = NULL, level = 0.95) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  st <- .standardize(X); Xs <- st$Xs; yc <- y - mean(y)     # ||x_j||^2 = n
  if (is.null(targets)) targets <- seq_len(p)

  # --- lasso principal, avec sigma estimé par scaled lasso (Sun-Zhang 2012) ---
  # On itère lambda proportionnel à sigma et sigma = ||residus||/sqrt(n-s) jusqu'à
  # stabilisation : ainsi lambda se cale sur le BRUIT, pas sur sd(y).
  sig <- if (is.null(sigma)) sd(yc) else sigma
  bs <- NULL
  for (iter in seq_len(20)) {
    lam <- if (is.null(lambda)) .lambda_theo(n, p, sig) else lambda
    fit <- lasso_fit(Xs, yc, lambda = lam, standardize = FALSE, intercept = FALSE, tol = 1e-10)
    bs <- fit$beta
    r <- yc - as.numeric(Xs %*% bs)
    s_hat <- sum(bs != 0)
    sig_new <- sqrt(sum(r^2) / max(n - s_hat, 1))
    if (!is.null(sigma) || !is.null(lambda) || abs(sig_new - sig) < 1e-3 * sig) { sig <- sig_new; break }
    sig <- sig_new
  }
  lambda <- if (is.null(lambda)) .lambda_theo(n, p, sig) else lambda
  if (is.null(lambda_node)) lambda_node <- .lambda_theo(n, p, 1)

  est <- se <- numeric(length(targets)); names(est) <- names(se) <- colnames(X)[targets]
  for (m in seq_along(targets)) {
    j <- targets[m]
    gam <- lasso_fit(Xs[, -j, drop = FALSE], Xs[, j], lambda = lambda_node,
                     standardize = FALSE, intercept = FALSE, tol = 1e-10)$beta
    tau <- Xs[, j] - as.numeric(Xs[, -j, drop = FALSE] %*% gam)   # score orthogonal
    denom <- sum(tau * Xs[, j])
    bd <- bs[j] + sum(tau * r) / denom
    se_s <- sig * sqrt(sum(tau^2)) / abs(denom)
    # retour à l'échelle d'origine (Xs_j = (x_j - center)/scale_j)
    est[m] <- bd / st$scale[j]; se[m] <- se_s / st$scale[j]
  }
  z <- qnorm(1 - (1 - level) / 2)
  list(estimate = est, se = se, lower = est - z * se, upper = est + z * se,
       beta_lasso = (bs / st$scale)[targets], sigma = sig, lambda = lambda)
}
