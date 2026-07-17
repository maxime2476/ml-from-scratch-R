# =============================================================================
# Module 46 — Analyse de survie (duree)
# Implemente les equations de derivations/46_survival.qmd. R base.
# La variable expliquee est un TEMPS jusqu'a un evenement, souvent CENSURE (on
# sait seulement qu'il depasse la duree observee). L'OLS/GLM ne gerent pas la
# censure ; Kaplan-Meier et Cox si.
# =============================================================================

#' Estimateur de survie de Kaplan-Meier
#'
#' \eqn{\hat S(t)=\prod_{t_i\le t}\bigl(1-d_i/n_i\bigr)}, ou \eqn{d_i} deces et
#' \eqn{n_i} sujets a risque en \eqn{t_i}. Gere la **censure** a droite (les
#' censures reduisent le risque sans compter comme evenement).
#'
#' @param time durees observees
#' @param event 1 = evenement, 0 = censure.
#' @return liste : `time` (temps d'evenement), `surv`, `n_risk`, `n_event`.
#' @export
kaplan_meier <- function(time, event) {
  ut <- sort(unique(time[event == 1])); S <- 1
  surv <- nr <- ne <- numeric(length(ut))
  for (i in seq_along(ut)) {
    t <- ut[i]; d <- sum(time == t & event == 1); risk <- sum(time >= t)
    S <- S * (1 - d / risk); surv[i] <- S; nr[i] <- risk; ne[i] <- d
  }
  list(time = ut, surv = surv, n_risk = nr, n_event = ne)
}

#' Test du log-rank (comparaison de courbes de survie)
#'
#' Compare les deces **observes** et **attendus** entre groupes a chaque temps
#' d'evenement ; \eqn{\chi^2} sous l'hypothese de survies egales.
#'
#' @param time,event durees et indicateurs
#' @param group facteur a 2 niveaux.
#' @return liste : `statistic`, `df`, `p_value`.
#' @export
logrank_test <- function(time, event, group) {
  g <- as.integer(as.factor(group)) - 1L; ut <- sort(unique(time[event == 1]))
  O1 <- E1 <- V <- 0
  for (t in ut) {
    nr <- sum(time >= t); d <- sum(time == t & event == 1)
    n1 <- sum(time >= t & g == 1); d1 <- sum(time == t & event == 1 & g == 1)
    O1 <- O1 + d1; E1 <- E1 + d * n1 / nr
    if (nr > 1) V <- V + d * (n1 / nr) * (1 - n1 / nr) * (nr - d) / (nr - 1)
  }
  stat <- (O1 - E1)^2 / V
  list(statistic = stat, df = 1, p_value = pchisq(stat, 1, lower.tail = FALSE))
}

#' Modele de Cox a risques proportionnels (vraisemblance partielle)
#'
#' Le risque instantane est \eqn{\lambda(t\mid x)=\lambda_0(t)\exp(x^\top\beta)}
#' (proportionnel : le risque de base \eqn{\lambda_0} n'a pas besoin d'etre
#' specifie). On estime \eqn{\beta} par la **vraisemblance partielle** de Cox
#' \eqn{\prod_{i:\text{evt}}\exp(x_i^\top\beta)/\sum_{j\in R(t_i)}\exp(x_j^\top\beta)}
#' (gestion des ex-aequo de Breslow), maximisee par Newton.
#'
#' @param time,event durees et indicateurs
#' @param X matrice de covariables.
#' @return liste : `coefficients`, `se`, `loglik`, `hazard_ratio`.
#' @export
cox_ph <- function(time, event, X) {
  X <- as.matrix(X); o <- order(time); time <- time[o]; event <- event[o]; X <- X[o, , drop = FALSE]
  ev <- which(event == 1); p <- ncol(X)
  neg_ll <- function(b) {
    eta <- as.numeric(X %*% b); -sum(sapply(ev, function(i) eta[i] - log(sum(exp(eta[time >= time[i]])))))
  }
  grad <- function(b) {
    eta <- as.numeric(X %*% b); g <- rep(0, p)
    for (i in ev) { risk <- time >= time[i]; w <- exp(eta[risk]); w <- w / sum(w)
      g <- g + X[i, ] - colSums(w * X[risk, , drop = FALSE]) }
    -g
  }
  opt <- optim(rep(0, p), neg_ll, grad, method = "BFGS", hessian = TRUE)
  vcov <- solve(opt$hessian); b <- opt$par; names(b) <- colnames(X)
  list(coefficients = b, se = sqrt(diag(vcov)), loglik = -opt$value, hazard_ratio = exp(b))
}
