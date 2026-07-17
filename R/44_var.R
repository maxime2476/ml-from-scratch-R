# =============================================================================
# Module 44 — Modeles vectoriels autoregressifs (VAR)
# Implemente les equations de derivations/44_var.qmd. R base.
# Le Module 31 etait UNIVARIE. Un VAR modelise PLUSIEURS series conjointement :
# chaque variable depend de ses retards ET de ceux des autres. Outils cles : la
# causalite de Granger et les fonctions de reponse impulsionnelle (IRF).
# =============================================================================

#' Estimation d'un VAR(p) (equation par equation, OLS)
#'
#' \eqn{Y_t=c+A_1Y_{t-1}+\dots+A_pY_{t-p}+\varepsilon_t}. Chaque equation est une
#' regression OLS sur la constante et les retards empiles.
#'
#' @param Y matrice n x k des series
#' @param p ordre.
#' @return objet `var_scratch` : `B` (coefficients), `Sigma`, `A` (liste des A_l), etc.
#' @export
var_fit <- function(Y, p = 1L) {
  Y <- as.matrix(Y); n <- nrow(Y); k <- ncol(Y)
  Z <- cbind(1, do.call(cbind, lapply(seq_len(p), function(l) Y[(p + 1 - l):(n - l), , drop = FALSE])))
  Yt <- Y[(p + 1):n, , drop = FALSE]
  B <- solve(crossprod(Z), crossprod(Z, Yt))                # (1 + k p) x k
  res <- Yt - Z %*% B; Sigma <- crossprod(res) / (nrow(Yt) - ncol(Z))
  A <- lapply(seq_len(p), function(l) t(B[(2 + (l - 1) * k):(1 + l * k), , drop = FALSE]))
  structure(list(B = B, Sigma = Sigma, A = A, c = B[1, ], p = p, k = k, res = res,
                 names = colnames(Y), Y = Y), class = "var_scratch")
}

#' Test de causalite de Granger
#'
#' \eqn{x} **cause au sens de Granger** \eqn{y} si ses retards ameliorent la
#' prediction de \eqn{y}. Test F comparant l'equation de \eqn{y} avec et sans les
#' retards de \eqn{x}.
#'
#' @param object objet `var_fit`
#' @param cause indice de la variable causante
#' @param effect indice de la variable expliquee.
#' @return liste : `statistic` (F), `df1`, `df2`, `p_value`.
#' @export
granger_test <- function(object, cause, effect) {
  Y <- object$Y; p <- object$p; k <- object$k; n <- nrow(Y)
  Z <- cbind(1, do.call(cbind, lapply(seq_len(p), function(l) Y[(p + 1 - l):(n - l), , drop = FALSE])))
  ye <- Y[(p + 1):n, effect]
  rss_f <- sum(lm.fit(Z, ye)$residuals^2)
  drop <- 1 + (0:(p - 1)) * k + cause                        # colonnes des retards de 'cause'
  rss_r <- sum(lm.fit(Z[, -drop, drop = FALSE], ye)$residuals^2)
  df1 <- p; df2 <- length(ye) - ncol(Z)
  Fst <- ((rss_r - rss_f) / df1) / (rss_f / df2)
  list(statistic = Fst, df1 = df1, df2 = df2, p_value = pf(Fst, df1, df2, lower.tail = FALSE))
}

#' Fonctions de reponse impulsionnelle (IRF) orthogonalisees
#'
#' Reponse dynamique de chaque variable a un choc unitaire (orthogonalise par
#' Cholesky de \eqn{\Sigma}) dans chaque variable, sur `h` periodes. Coeur de
#' l'analyse macro-econometrique (propagation des chocs).
#'
#' @param object objet `var_fit`
#' @param h horizon.
#' @return tableau (h+1) x k x k : `irf[t, reponse, choc]`.
#' @export
var_irf <- function(object, h = 10L) {
  k <- object$k; p <- object$p; A <- object$A
  Phi <- vector("list", h + 1); Phi[[1]] <- diag(k)         # MA : Phi_0 = I
  for (i in seq_len(h)) {
    S <- matrix(0, k, k)
    for (l in seq_len(min(i, p))) S <- S + A[[l]] %*% Phi[[i - l + 1]]
    Phi[[i + 1]] <- S
  }
  P <- t(chol(object$Sigma))                                 # Cholesky inferieur
  out <- array(0, c(h + 1, k, k))
  for (i in seq_len(h + 1)) out[i, , ] <- Phi[[i]] %*% P
  out
}
