# =============================================================================
# Module 20 — Régression quantile
# Implémente les équations de derivations/20_quantile.qmd. R base + Module 0.
# =============================================================================

#' Perte pinball (fonction « check », éq. 20.1)
#'
#' \eqn{\rho_\tau(u)=u(\tau-\mathbf 1\{u<0\})}. Vectorisée.
#'
#' @param u résidu(s).
#' @param tau niveau de quantile dans (0,1).
#' @return la (les) perte(s) pinball.
#' @export
pinball_loss <- function(u, tau) u * (tau - (u < 0))

#' Régression quantile par IRLS (éq. 20.2)
#'
#' Minimise la perte pinball par moindres carrés pondérés itérés (majoration de
#' Hunter-Lange). Réutilise la QR pondérée du Module 0. Reproche `quantreg::rq`.
#'
#' @param formula formule façon `lm`.
#' @param data data.frame.
#' @param tau niveau de quantile (défaut 0.5 = médiane / LAD).
#' @param maxit itérations maximales.
#' @param tol tolérance d'arrêt (variation des coefficients).
#' @param eps perturbation de la majoration (évite la division par 0).
#' @return liste : `coefficients`, `tau`, `fitted`, `residuals`, `loss`, `iter`.
#' @export
qreg_fit <- function(formula, data, tau = 0.5, maxit = 200L, tol = 1e-8, eps = 1e-6) {
  mf <- model.frame(formula, data)
  y <- as.numeric(model.response(mf)); X <- model.matrix(attr(mf, "terms"), mf)
  n <- nrow(X); p <- ncol(X)
  beta <- solve_ls_qr(X, y)$coefficients                 # départ OLS
  for (it in seq_len(maxit)) {
    r <- as.numeric(y - X %*% beta)
    a <- ifelse(r >= 0, tau, 1 - tau)
    w <- a / pmax(abs(r), eps)                            # poids IRLS (éq. 20.2)
    sw <- sqrt(w)
    beta_new <- solve_ls_qr(X * sw, y * sw)$coefficients
    if (max(abs(beta_new - beta)) < tol) { beta <- beta_new; break }
    beta <- beta_new
  }
  names(beta) <- colnames(X)
  fitted <- as.numeric(X %*% beta); resid <- y - fitted
  list(coefficients = beta, tau = tau, fitted = fitted, residuals = resid,
       loss = sum(pinball_loss(resid, tau)), iter = it)
}
