# =============================================================================
# Module 45 — Volatilite conditionnelle (GARCH)
# Implemente les equations de derivations/45_garch.qmd. R base.
# Les rendements financiers ont une variance qui VARIE dans le temps : les
# periodes agitees s'enchainent (regroupement de volatilite). Le GARCH modelise
# cette variance conditionnelle comme fonction du passe des chocs et de la variance.
# =============================================================================

#' Test ARCH-LM (heteroscedasticite conditionnelle)
#'
#' Regresse \eqn{\hat\varepsilon_t^2} sur ses \eqn{q} retards ; \eqn{nR^2\sim
#' \chi^2_q} sous absence d'effet ARCH. Rejet = la volatilite se regroupe.
#'
#' @param x serie (rendements) ; @param q nombre de retards.
#' @return liste : `statistic`, `df`, `p_value`.
#' @export
arch_lm_test <- function(x, q = 5L) {
  e2 <- (x - mean(x))^2; n <- length(e2)
  Z <- cbind(1, sapply(seq_len(q), function(l) c(rep(0, l), e2[seq_len(n - l)])))
  r <- lm.fit(Z, e2)$residuals; r2 <- 1 - sum(r^2) / sum((e2 - mean(e2))^2)
  stat <- (n - q) * r2
  list(statistic = stat, df = q, p_value = pchisq(stat, q, lower.tail = FALSE))
}

#' Estimation d'un GARCH(1,1) par (quasi-)maximum de vraisemblance
#'
#' \eqn{x_t=\sigma_t z_t}, \eqn{z_t\sim(0,1)}, avec la variance conditionnelle
#' \eqn{\sigma_t^2=\omega+\alpha x_{t-1}^2+\beta\sigma_{t-1}^2}. On maximise la
#' log-vraisemblance gaussienne (paramétrage assurant \eqn{\omega>0}, \eqn{\alpha,
#' \beta\ge0}, \eqn{\alpha+\beta<1} : stationnarite).
#'
#' @param x serie (rendements centres) ; @param maxit iterations de l'optimiseur.
#' @return liste : `omega`, `alpha`, `beta`, `sigma` (volatilite conditionnelle), `loglik`.
#' @export
garch_fit <- function(x, maxit = 500L) {
  x <- x - mean(x); n <- length(x); v <- var(x)
  cond_var <- function(w, a, b) {
    s2 <- numeric(n); s2[1] <- v
    for (t in 2:n) s2[t] <- w + a * x[t - 1]^2 + b * s2[t - 1]
    s2
  }
  nll <- function(par) {
    w <- exp(par[1]); a <- plogis(par[2]) * 0.999; b <- plogis(par[3]) * (0.999 - a)
    s2 <- cond_var(w, a, b); 0.5 * sum(log(s2) + x^2 / s2)
  }
  o <- optim(c(log(0.05), qlogis(0.1), qlogis(0.5)), nll, method = "BFGS", control = list(maxit = maxit))
  w <- exp(o$par[1]); a <- plogis(o$par[2]) * 0.999; b <- plogis(o$par[3]) * (0.999 - a)
  list(omega = w, alpha = a, beta = b, sigma = sqrt(cond_var(w, a, b)),
       persistence = a + b, loglik = -o$value)
}
