# =============================================================================
# Module 30 — Variables dependantes limitees
# Implemente les equations de derivations/30_limited_dependent.qmd. R base.
# Quand la variable expliquee est binaire (probit), censuree (Tobit) ou observee
# seulement pour un sous-echantillon selectionne (Heckman), l'OLS est biaise :
# il faut modeliser la structure de la reponse par maximum de vraisemblance.
# =============================================================================

.mm <- function(formula, data) {
  mf <- model.frame(formula, data, na.action = na.pass); tt <- attr(mf, "terms")
  list(y = as.numeric(model.response(mf)), X = model.matrix(tt, mf))
}

#' Probit (reponse binaire) par IRLS (lien probit)
#'
#' \eqn{P(y=1)=\Phi(X\beta)}. Maximum de vraisemblance par moindres carres
#' ponderes iteres (IRLS, Module 3) avec lien probit : poids
#' \eqn{w=\phi^2/[\Phi(1-\Phi)]}, reponse de travail \eqn{z=\eta+(y-\Phi)/\phi}.
#'
#' @param formula formule 
#' @param data data.frame 
#' @param maxit iterations max
#' @return liste : `coefficients`, `vcov`, `se`, `fitted`, `loglik`.
#' @export
probit <- function(formula, data, maxit = 50L) {
  d <- .mm(formula, data); X <- d$X; y <- d$y; b <- rep(0, ncol(X))
  for (it in seq_len(maxit)) {
    eta <- as.numeric(X %*% b); p <- pmin(pmax(pnorm(eta), 1e-8), 1 - 1e-8); ph <- dnorm(eta)
    w <- ph^2 / (p * (1 - p)); z <- eta + (y - p) / ph
    bn <- solve(crossprod(X * w, X), crossprod(X * w, z))
    if (max(abs(bn - b)) < 1e-10) { b <- bn; break }; b <- bn
  }
  eta <- as.numeric(X %*% b); p <- pmin(pmax(pnorm(eta), 1e-8), 1 - 1e-8); ph <- dnorm(eta)
  w <- ph^2 / (p * (1 - p)); vcov <- solve(crossprod(X * w, X))
  beta <- as.numeric(b); names(beta) <- colnames(X); dimnames(vcov) <- list(colnames(X), colnames(X))
  list(coefficients = beta, vcov = vcov, se = sqrt(diag(vcov)),
       fitted = p, loglik = sum(y * log(p) + (1 - y) * log(1 - p)))
}

#' Regression Tobit (reponse censuree) par maximum de vraisemblance
#'
#' Modele a variable latente \eqn{y^*=X\beta+\varepsilon}, \eqn{\varepsilon\sim
#' \mathcal N(0,\sigma^2)}, observe \eqn{y=\max(L,y^*)}. Vraisemblance : densite
#' normale pour les non censures, \eqn{\Phi((L-X\beta)/\sigma)} pour les censures.
#'
#' @param formula formule 
#' @param data data.frame 
#' @param left seuil de censure
#' @return liste : `coefficients`, `sigma`, `se`, `loglik`.
#' @export
tobit_fit <- function(formula, data, left = 0) {
  d <- .mm(formula, data); X <- d$X; y <- d$y; k <- ncol(X)
  nll <- function(par) {
    b <- par[1:k]; s <- exp(par[k + 1]); mu <- as.numeric(X %*% b); unc <- y > left
    -(sum(dnorm((y[unc] - mu[unc]) / s, log = TRUE) - log(s)) +
      sum(pnorm((left - mu[!unc]) / s, log.p = TRUE)))
  }
  init <- c(coef(lm(y ~ X - 1)), log(sd(y)))
  opt <- optim(init, nll, method = "BFGS", hessian = TRUE)
  vcov <- solve(opt$hessian); se <- sqrt(diag(vcov))
  beta <- opt$par[1:k]; names(beta) <- colnames(X)
  sigma <- exp(opt$par[k + 1])
  list(coefficients = beta, sigma = sigma, se = se[1:k],
       se_logsigma = se[k + 1], loglik = -opt$value)
}

#' Modele de selection de Heckman (estimation en deux etapes)
#'
#' Corrige le **biais de selection** : \eqn{y} n'est observe que si \eqn{d=1}
#' (equation de selection). Etape 1 : probit de \eqn{d} sur \eqn{Z_s}, ratio de
#' Mills inverse \eqn{\hat\lambda=\phi(Z_s\hat\gamma)/\Phi(Z_s\hat\gamma)}.
#' Etape 2 : OLS de \eqn{y} sur \eqn{[X_o,\hat\lambda]} sur les observations
#' selectionnees. Le coefficient de \eqn{\hat\lambda} vaut \eqn{\rho\sigma}
#' (nul \eqn{\iff} pas de biais de selection).
#'
#' @param selection formule de selection (LHS binaire).
#' @param outcome formule de resultat (LHS avec NA hors selection).
#' @param data data.frame.
#' @return liste : `gamma` (selection), `beta` (resultat, dont `imr`), `se`.
#' @export
heckman <- function(selection, outcome, data) {
  ds <- .mm(selection, data); Zs <- ds$X; dsel <- ds$y
  pr <- probit(selection, data)
  eta <- as.numeric(Zs %*% pr$coefficients); imr <- dnorm(eta) / pnorm(eta)
  do <- .mm(outcome, data); Xo <- do$X; yo <- do$y
  obs <- dsel == 1 & !is.na(yo)
  Xa <- cbind(Xo[obs, , drop = FALSE], imr = imr[obs])
  beta <- solve(crossprod(Xa), crossprod(Xa, yo[obs]))
  r <- yo[obs] - Xa %*% beta; s2 <- sum(r^2) / (sum(obs) - ncol(Xa))
  vcov <- s2 * solve(crossprod(Xa))                       # SE naives (voir remarque)
  b <- as.numeric(beta); names(b) <- c(colnames(Xo), "imr")
  list(gamma = pr$coefficients, beta = b, se = sqrt(diag(vcov)),
       lambda = imr, n_obs = sum(obs))
}
