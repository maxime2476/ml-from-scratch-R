# =============================================================================
# Module 29 — Tests de specification et diagnostics
# Implemente les equations de derivations/29_diagnostics.qmd. R base.
# Le pendant "detection" des corrections du Module 2 (sandwich) : d'abord
# DIAGNOSTIQUER heteroscedasticite, autocorrelation, endogeneite, mauvaise
# specification, puis y REMEDIER (FGLS ici, IV au Module 5, HAC au Module 2).
# =============================================================================

# utilitaire interne : ajuste l'OLS et renvoie design, residus, etc.
.ols_parts <- function(formula, data) {
  mf <- model.frame(formula, data); tt <- attr(mf, "terms")
  y <- as.numeric(model.response(mf)); X <- model.matrix(tt, mf)
  n <- nrow(X); p <- ncol(X)
  beta <- as.numeric(solve(crossprod(X), crossprod(X, y)))
  fitted <- as.numeric(X %*% beta); e <- y - fitted
  list(y = y, X = X, beta = beta, fitted = fitted, e = e, n = n, p = p, terms = tt, mf = mf)
}
.aux_nr2 <- function(response, Z) {           # n*R^2 d'une regression auxiliaire
  Z <- cbind(1, Z); n <- length(response)
  b <- solve(crossprod(Z), crossprod(Z, response)); r <- response - Z %*% b
  r2 <- 1 - sum(r^2) / sum((response - mean(response))^2)
  list(stat = n * r2, r2 = r2)
}

#' Test de Breusch-Pagan (heteroscedasticite)
#'
#' Version studentisee de Koenker : \eqn{n R^2} de la regression des residus au
#' carre sur les regresseurs. \eqn{H_0} : homoscedasticite. Statistique
#' \eqn{\sim \chi^2_{p-1}}.
#'
#' @param formula formule du modele
#' @param data data.frame
#' @return liste : `statistic`, `df`, `p_value`.
#' @export
bp_test <- function(formula, data) {
  o <- .ols_parts(formula, data)
  a <- .aux_nr2(o$e^2, o$X[, -1, drop = FALSE]); df <- o$p - 1
  list(statistic = a$stat, df = df, p_value = pchisq(a$stat, df, lower.tail = FALSE))
}

#' Test de White (heteroscedasticite generale)
#'
#' Breusch-Pagan augmente des carres et produits croises des regresseurs :
#' detecte toute forme d'heteroscedasticite liee a \eqn{X}.
#'
#' @param formula,data cf. `bp_test`.
#' @return liste : `statistic`, `df`, `p_value`.
#' @export
white_test <- function(formula, data) {
  o <- .ols_parts(formula, data); Xr <- o$X[, -1, drop = FALSE]; k <- ncol(Xr)
  terms_ <- Xr
  for (i in seq_len(k)) for (j in i:k) terms_ <- cbind(terms_, Xr[, i] * Xr[, j])  # carres + croises
  terms_ <- terms_[, !duplicated(round(t(terms_), 10)), drop = FALSE]
  a <- .aux_nr2(o$e^2, terms_); df <- ncol(terms_)
  list(statistic = a$stat, df = df, p_value = pchisq(a$stat, df, lower.tail = FALSE))
}

#' Statistique de Durbin-Watson (autocorrelation d'ordre 1)
#'
#' \eqn{DW = \sum_{t\ge2}(e_t-e_{t-1})^2/\sum_t e_t^2 \approx 2(1-\hat\rho)}.
#' Proche de 2 : pas d'autocorrelation ; proche de 0 : positive.
#'
#' @param formula,data cf. `bp_test`.
#' @return liste : `statistic`, `rho` (autocorrelation implicite).
#' @export
dw_test <- function(formula, data) {
  e <- .ols_parts(formula, data)$e
  dw <- sum(diff(e)^2) / sum(e^2)
  list(statistic = dw, rho = 1 - dw / 2)
}

#' Test de Breusch-Godfrey (autocorrelation d'ordre p, LM)
#'
#' Regression des residus sur les regresseurs ET \eqn{p} retards des residus ;
#' \eqn{n R^2 \sim \chi^2_p}. Plus general que Durbin-Watson (ordre eleve,
#' regresseurs retardes admis).
#'
#' @param formula,data cf. `bp_test` 
#' @param order ordre du retard
#' @return liste : `statistic`, `df`, `p_value`.
#' @export
bg_test <- function(formula, data, order = 1) {
  o <- .ols_parts(formula, data); e <- o$e; n <- o$n
  L <- sapply(seq_len(order), function(k) { v <- c(rep(0, k), e[seq_len(n - k)]); v })
  a <- .aux_nr2(e, cbind(o$X[, -1, drop = FALSE], L))
  # n*R^2 sur la regression complete e ~ X + retards ; df = order
  Z <- cbind(o$X, L); b <- solve(crossprod(Z), crossprod(Z, e)); r <- e - Z %*% b
  r2 <- 1 - sum(r^2) / sum(e^2); stat <- n * r2
  list(statistic = stat, df = order, p_value = pchisq(stat, order, lower.tail = FALSE))
}

#' Test RESET de Ramsey (forme fonctionnelle)
#'
#' Ajoute des puissances des valeurs ajustees \eqn{\hat y^2,\hat y^3,\dots} au
#' modele et teste (F) leur significativite jointe. Rejet = non-linearite
#' negligee.
#'
#' @param formula,data cf. `bp_test` 
#' @param powers puissances de \eqn{\hat y}
#' @return liste : `statistic` (F), `df1`, `df2`, `p_value`.
#' @export
reset_test <- function(formula, data, powers = 2:3) {
  o <- .ols_parts(formula, data)
  P <- sapply(powers, function(k) o$fitted^k)
  Za <- cbind(o$X, P); ba <- solve(crossprod(Za), crossprod(Za, o$y))
  rss_a <- sum((o$y - Za %*% ba)^2); rss_0 <- sum(o$e^2)
  df1 <- length(powers); df2 <- o$n - ncol(Za)
  Fst <- ((rss_0 - rss_a) / df1) / (rss_a / df2)
  list(statistic = Fst, df1 = df1, df2 = df2, p_value = pf(Fst, df1, df2, lower.tail = FALSE))
}

#' Test de normalite de Jarque-Bera
#'
#' \eqn{JB = \tfrac n6\bigl(S^2 + \tfrac14(K-3)^2\bigr)}, \eqn{S} asymetrie,
#' \eqn{K} aplatissement. \eqn{\sim\chi^2_2} sous normalite.
#'
#' @param x vecteur (typiquement des residus).
#' @return liste : `statistic`, `skewness`, `kurtosis`, `p_value`.
#' @export
jarque_bera <- function(x) {
  x <- x - mean(x); n <- length(x); m2 <- mean(x^2)
  S <- mean(x^3) / m2^1.5; K <- mean(x^4) / m2^2
  jb <- n / 6 * (S^2 + (K - 3)^2 / 4)
  list(statistic = jb, skewness = S, kurtosis = K, p_value = pchisq(jb, 2, lower.tail = FALSE))
}

#' Test d'endogeneite de Durbin-Wu-Hausman (regression augmentee)
#'
#' Regresse les regresseurs suspects sur les instruments, recupere les residus,
#' les ajoute a l'OLS et teste leur significativite jointe (F). Rejet = OLS
#' biaise, l'IV (Module 5) est requis.
#'
#' @param y reponse 
#' @param X design (constante + regresseurs, endogenes inclus)
#' @param Z instruments (constante + exogenes + exclus).
#' @param endog indices des colonnes ENDOGENES de X.
#' @return liste : `statistic` (F), `df1`, `df2`, `p_value`.
#' @export
dwh_test <- function(y, X, Z, endog) {
  X <- as.matrix(X); Z <- as.matrix(Z); y <- as.numeric(y); n <- nrow(X)
  V <- sapply(endog, function(j) X[, j] - Z %*% solve(crossprod(Z), crossprod(Z, X[, j])))
  Xa <- cbind(X, V); ba <- solve(crossprod(Xa), crossprod(Xa, y))
  rss_a <- sum((y - Xa %*% ba)^2)
  b0 <- solve(crossprod(X), crossprod(X, y)); rss_0 <- sum((y - X %*% b0)^2)
  df1 <- length(endog); df2 <- n - ncol(Xa)
  Fst <- ((rss_0 - rss_a) / df1) / (rss_a / df2)
  list(statistic = Fst, df1 = df1, df2 = df2, p_value = pf(Fst, df1, df2, lower.tail = FALSE))
}

#' Test de suridentification de Sargan
#'
#' \eqn{n R^2} de la regression des residus 2SLS sur TOUS les instruments ;
#' \eqn{\sim\chi^2_{q-k}} (q instruments, k regresseurs). Rejet = instruments
#' invalides (correles au terme d'erreur).
#'
#' @param y reponse 
#' @param X design (endogenes inclus) 
#' @param Z instruments
#' @return liste : `statistic`, `df`, `p_value`.
#' @export
sargan_test <- function(y, X, Z) {
  X <- as.matrix(X); Z <- as.matrix(Z); y <- as.numeric(y); n <- nrow(X)
  Xhat <- Z %*% solve(crossprod(Z), crossprod(Z, X))            # projection sur Z
  b2sls <- solve(crossprod(Xhat), crossprod(Xhat, y))
  r <- as.numeric(y - X %*% b2sls)
  r2 <- 1 - sum((r - Z %*% solve(crossprod(Z), crossprod(Z, r)))^2) / sum(r^2)
  stat <- n * r2; df <- ncol(Z) - ncol(X)
  list(statistic = stat, df = df, p_value = pchisq(stat, df, lower.tail = FALSE))
}

#' Moindres carres generalises FAISABLES (FGLS) pour heteroscedasticite
#'
#' Modele de variance : \eqn{\log\hat\varepsilon^2 = X\gamma}. On en deduit des
#' poids \eqn{\hat w_i = 1/\exp(X_i\hat\gamma)} et l'on refait une WLS (Module 2).
#' Plus efficace que l'OLS sous heteroscedasticite bien specifiee.
#'
#' @param formula,data cf. `bp_test`.
#' @return liste : `coefficients`, `weights`, `se` (classiques WLS).
#' @export
fgls <- function(formula, data) {
  o <- .ols_parts(formula, data)
  g <- solve(crossprod(o$X), crossprod(o$X, log(o$e^2 + 1e-12)))  # modele de variance
  w <- as.numeric(1 / exp(o$X %*% g))
  sw <- sqrt(w); Xw <- o$X * sw; yw <- o$y * sw
  b <- solve(crossprod(Xw), crossprod(Xw, yw))
  r <- yw - Xw %*% b; s2 <- sum(r^2) / (o$n - o$p)
  vcov <- s2 * solve(crossprod(Xw))
  beta <- as.numeric(b); names(beta) <- colnames(o$X)
  list(coefficients = beta, weights = w, se = sqrt(diag(vcov)))
}
