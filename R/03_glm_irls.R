# =============================================================================
# Module 3 — GLM par moindres carrés pondérés itérés (IRLS)
# Implémente les équations de derivations/03_glm_irls.qmd.
# Familles : binomiale (logit) et Poisson (log), liens canoniques.
# Chaque itération IRLS résout une WLS via la QR du Module 0 (solve_ls_qr).
# =============================================================================

# ---- Familles (lien canonique) ----------------------------------------------
# Chaque famille fournit : linkinv (mu <- eta), mu_eta (dmu/deta), variance
# V(mu), dev_resids (contributions à la déviance, éq. 3.14), loglik, mustart.
.ylogy <- function(y, mu) ifelse(y == 0, 0, y * log(y / mu))

.family_binomial <- list(
  name = "binomial",
  linkfun  = function(mu)  log(mu / (1 - mu)),
  linkinv  = function(eta) plogis(eta),                 # sigmoïde
  mu_eta   = function(mu)  mu * (1 - mu),                # dmu/deta = W (logit)
  variance = function(mu)  mu * (1 - mu),
  dev      = function(y, mu) 2 * (.ylogy(y, mu) + .ylogy(1 - y, 1 - mu)),
  loglik   = function(y, mu) sum(y * log(mu) + (1 - y) * log(1 - mu)),
  mustart  = function(y) (y + 0.5) / 2
)

.family_poisson <- list(
  name = "poisson",
  linkfun  = function(mu)  log(mu),
  linkinv  = function(eta) exp(eta),
  mu_eta   = function(mu)  mu,                            # dmu/deta = mu (log)
  variance = function(mu)  mu,
  dev      = function(y, mu) 2 * (.ylogy(y, mu) - (y - mu)),
  loglik   = function(y, mu) sum(dpois(y, mu, log = TRUE)),
  mustart  = function(y) y + 0.1
)

.get_family <- function(family) {
  switch(family,
    binomial = .family_binomial,
    poisson  = .family_poisson,
    stop("family doit être 'binomial' ou 'poisson'."))
}

#' Ajustement d'un GLM par IRLS
#'
#' Boucle IRLS (éq. 3.7-3.9) : à chaque étape, poids \eqn{W}, réponse de travail
#' \eqn{z} (éq. 3.8), puis WLS de z sur X (QR du Module 0). Convergence sur la
#' variation relative de la déviance (critère de `glm`). Variance par
#' l'information de Fisher \eqn{(X^TWX)^{-1}} (éq. 3.10).
#'
#' @param formula formule façon `glm`.
#' @param data data.frame.
#' @param family "binomial" (logit) ou "poisson" (log).
#' @param maxit itérations maximales (défaut 25, comme `glm`).
#' @param epsilon tolérance de convergence (défaut 1e-8, comme `glm`).
#' @return objet `glm_irls` : `coefficients`, `vcov`, `se`, `fitted`,
#'   `linear.predictors`, `deviance`, `null.deviance`, `loglik`, `iter`,
#'   `df.residual`, `rank`, `family`, `weights`, `model_matrix`, `response`.
glm_irls <- function(formula, data, family = c("binomial", "poisson"),
                     maxit = 25L, epsilon = 1e-8) {
  family <- match.arg(family)
  fam <- .get_family(family)
  mf <- model.frame(formula, data)
  tt <- attr(mf, "terms")
  y  <- as.numeric(model.response(mf))
  X  <- model.matrix(tt, mf)
  n <- nrow(X); p <- ncol(X)

  mu  <- fam$mustart(y)
  eta <- fam$linkfun(mu)
  devold <- sum(fam$dev(y, mu))
  beta <- rep(0, p)
  W_last <- NULL; R1 <- NULL
  conv <- FALSE

  for (iter in seq_len(maxit)) {
    W_diag  <- fam$mu_eta(mu)^2 / fam$variance(mu)       # poids IRLS (éq. 3.9)
    z <- eta + (y - mu) / fam$mu_eta(mu)                  # réponse de travail (3.8)
    sw <- sqrt(W_diag)
    fit <- solve_ls_qr(X * sw, z * sw)                    # WLS via QR (Module 0)
    beta <- fit$coefficients
    W_last <- W_diag; R1 <- fit$R                         # info de Fisher de cette étape
    eta <- as.numeric(X %*% beta)
    mu  <- fam$linkinv(eta)
    dev <- sum(fam$dev(y, mu))
    if (abs(dev - devold) / (abs(dev) + 0.1) < epsilon) { conv <- TRUE; break }
    devold <- dev
  }
  if (!conv) warning("IRLS n'a pas convergé en ", maxit, " itérations.")

  # Information de Fisher et variance (éq. 3.10), depuis le facteur R de la
  # dernière itération IRLS (convention de `glm` : poids au mu pré-mise-à-jour).
  Rinv <- backsolve(R1, diag(p))
  vcov <- Rinv %*% t(Rinv)                                # (X^T W X)^{-1}
  dimnames(vcov) <- list(colnames(X), colnames(X))
  names(beta) <- colnames(X)

  has_int <- attr(tt, "intercept") == 1
  mu_null <- if (has_int) fam$linkinv(fam$linkfun(mean(y))) else fam$linkinv(0)
  null.deviance <- if (has_int) sum(fam$dev(y, rep(mean(y), n))) else NA_real_

  structure(list(
    coefficients = beta, vcov = vcov, se = sqrt(diag(vcov)),
    fitted = mu, linear.predictors = eta,
    deviance = dev, null.deviance = null.deviance,
    loglik = fam$loglik(y, mu), iter = iter, df.residual = n - p, rank = p, n = n,
    family = family, weights = W_last, converged = conv,
    model_matrix = X, response = y, terms = tt
  ), class = "glm_irls")
}

#' Test de Wald (éq. 3.11) pour H0 : R beta = r
#'
#' @param fit objet `glm_irls`.
#' @param R matrice q x p des restrictions.
#' @param r vecteur cible (défaut zéro).
#' @return liste : `statistic`, `df`, `p_value`.
wald_test <- function(fit, R, r = rep(0, nrow(R))) {
  R <- as.matrix(R)
  Rb_r <- as.numeric(R %*% fit$coefficients - r)
  mid  <- R %*% fit$vcov %*% t(R)
  stat <- as.numeric(t(Rb_r) %*% solve(mid, Rb_r))
  q <- nrow(R)
  list(statistic = stat, df = q, p_value = pchisq(stat, q, lower.tail = FALSE))
}

#' Test du rapport de vraisemblance (éq. 3.12) pour modèles emboîtés
#'
#' \eqn{\mathrm{LR} = D_{\text{réduit}} - D_{\text{complet}}}.
#'
#' @param fit_full objet `glm_irls` du modèle complet.
#' @param fit_reduced objet `glm_irls` du modèle réduit (emboîté).
#' @return liste : `statistic`, `df`, `p_value`.
lr_test <- function(fit_full, fit_reduced) {
  stat <- fit_reduced$deviance - fit_full$deviance
  df <- fit_full$rank - fit_reduced$rank
  list(statistic = stat, df = df, p_value = pchisq(stat, df, lower.tail = FALSE))
}

#' Test du score / Rao (éq. 3.13) pour modèles emboîtés
#'
#' Statistique de Rao = réduction de la somme des carrés pondérée obtenue en
#' régressant les résidus de travail du modèle réduit sur le design complet, avec
#' les poids de travail du réduit (une itération de score de Fisher). C'est la
#' forme équivalente à \eqn{U(\tilde\beta)^T \mathcal I(\tilde\beta)^{-1}
#' U(\tilde\beta)} quand le réduit satisfait ses équations normales, et c'est
#' l'algorithme exact de `anova.glm(test = "Rao")`.
#'
#' @param fit_full objet `glm_irls` du modèle complet (fournit le design X).
#' @param fit_reduced objet `glm_irls` du modèle réduit (emboîté).
#' @return liste : `statistic`, `df`, `p_value`.
score_test <- function(fit_full, fit_reduced) {
  fam <- .get_family(fit_full$family)
  X <- fit_full$model_matrix
  y <- fit_full$response
  mu_r <- fit_reduced$fitted
  r <- (y - mu_r) / fam$mu_eta(mu_r)          # résidus de travail du réduit
  w <- fit_reduced$weights                     # poids de travail du réduit

  # Régression pondérée de r sur le design complet (glm gaussien = WLS).
  sw <- sqrt(w)
  bz <- solve_ls_qr(X * sw, r * sw)$coefficients
  fitted_r <- as.numeric(X %*% bz)
  wdev <- sum(w * (r - fitted_r)^2)            # déviance pondérée
  icpt <- attr(fit_full$terms, "intercept") == 1
  wmean <- if (icpt) sum(w * r) / sum(w) else 0
  wnull <- sum(w * (r - wmean)^2)              # déviance nulle pondérée
  stat <- wnull - wdev                          # Rao = null.deviance - deviance
  df <- fit_full$rank - fit_reduced$rank
  list(statistic = stat, df = df, p_value = pchisq(stat, df, lower.tail = FALSE))
}
