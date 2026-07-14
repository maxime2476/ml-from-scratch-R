# =============================================================================
# Module 33 — Methodes de Monte Carlo par chaines de Markov (MCMC)
# Implemente les equations de derivations/33_mcmc.qmd. R base.
# La voie COMPUTATIONNELLE de l'inference bayesienne : quand le posterieur n'a
# pas de forme fermee (contrairement aux Modules 14/27), on l'ECHANTILLONNE en
# construisant une chaine de Markov dont la loi stationnaire EST le posterieur.
# =============================================================================

#' Echantillonneur de Metropolis-Hastings (marche aleatoire)
#'
#' Propose \eqn{x'=x+\mathcal N(0,\text{psd}^2)} et l'accepte avec probabilite
#' \eqn{\min(1,\pi(x')/\pi(x))}. La chaine converge vers la loi cible \eqn{\pi}
#' (connue a une constante pres). Fonctionne pour toute cible via sa log-densite.
#'
#' @param log_target fonction `x -> log pi(x)` (a une constante additive pres).
#' @param init point de depart (scalaire ou vecteur).
#' @param proposal_sd ecart-type de la proposition (par composante).
#' @param n_iter nombre d'iterations.
#' @return liste : `chain` (n_iter x d), `accept_rate`.
#' @export
metropolis_hastings <- function(log_target, init, proposal_sd, n_iter = 10000L) {
  x <- init; k <- length(init); chain <- matrix(0, n_iter, k); acc <- 0
  lp <- log_target(x)
  for (i in seq_len(n_iter)) {
    prop <- x + rnorm(k, 0, proposal_sd); lpp <- log_target(prop)
    if (is.finite(lpp) && log(runif(1)) < lpp - lp) { x <- prop; lp <- lpp; acc <- acc + 1 }
    chain[i, ] <- x
  }
  list(chain = chain, accept_rate = acc / n_iter)
}

#' Gibbs pour la regression lineaire bayesienne (prior non informatif)
#'
#' Alterne les tirages des lois CONDITIONNELLES conjuguees :
#' \eqn{\beta\mid\sigma^2,y\sim\mathcal N(\hat\beta_{OLS},\sigma^2(X^\top X)^{-1})}
#' et \eqn{\sigma^2\mid\beta,y\sim\text{Inv-Gamma}(n/2,\;\|y-X\beta\|^2/2)}. La loi
#' stationnaire est le posterieur joint ; sa moyenne coincide avec l'OLS (M1).
#'
#' @param X design 
#' @param y reponse 
#' @param n_iter iterations 
#' @param burn rodage
#' @return liste : `beta` (echantillons apres burn), `sigma2`.
#' @export
gibbs_linreg <- function(X, y, n_iter = 5000L, burn = 1000L) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  XtXi <- solve(crossprod(X)); bhat <- as.numeric(XtXi %*% crossprod(X, y)); L <- t(chol(XtXi))
  s2 <- 1; B <- matrix(0, n_iter, p); S <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    b <- bhat + sqrt(s2) * as.numeric(L %*% rnorm(p))
    r <- y - X %*% b; s2 <- 1 / rgamma(1, n / 2, sum(r^2) / 2)
    B[i, ] <- b; S[i] <- s2
  }
  keep <- (burn + 1):n_iter
  list(beta = B[keep, , drop = FALSE], sigma2 = S[keep])
}

#' Diagnostic de convergence de Gelman-Rubin (statistique R-hat)
#'
#' Compare la variance INTER-chaines a la variance INTRA-chaine :
#' \eqn{\hat R=\sqrt{\hat V/W}}, \eqn{\hat V=\frac{n-1}n W+\frac Bn}. \eqn{\hat R\to
#' 1} a convergence ; \eqn{>1.1} signale une non-convergence.
#'
#' @param chains liste de chaines (vecteurs de meme longueur).
#' @return la statistique R-hat.
#' @export
gelman_rubin <- function(chains) {
  m <- length(chains); n <- length(chains[[1]])
  means <- vapply(chains, mean, numeric(1)); vars <- vapply(chains, var, numeric(1))
  B <- n * var(means); W <- mean(vars)
  Vhat <- (n - 1) / n * W + (1 + 1 / m) * B / n
  sqrt(Vhat / W)
}

#' Taille d'echantillon effective (ESS)
#'
#' \eqn{\mathrm{ESS}=n/(1+2\sum_k\hat\rho_k)}, somme des autocorrelations tronquee
#' au premier terme negatif (sequence positive initiale de Geyer). Mesure combien
#' de tirages INDEPENDANTS equivalent a la chaine correlee.
#'
#' @param x chaine (vecteur).
#' @return la taille d'echantillon effective.
#' @export
ess <- function(x) {
  n <- length(x); x <- x - mean(x)
  rho <- as.numeric(acf(x, lag.max = min(n - 1, 1000), plot = FALSE)$acf)[-1]
  s <- 0; for (k in seq_along(rho)) { if (rho[k] < 0) break; s <- s + rho[k] }
  n / (1 + 2 * s)
}
