# =============================================================================
# Module 31 — Series temporelles
# Implemente les equations de derivations/31_time_series.qmd. R base.
# La dependance temporelle invalide l'hypothese i.i.d. : on modelise la
# structure d'autocorrelation (AR/MA/ARMA), on la teste (Ljung-Box), et on
# distingue serie stationnaire d'une racine unitaire (Dickey-Fuller augmente).
# =============================================================================

#' Fonction d'autocorrelation (ACF)
#'
#' \eqn{\hat\rho_k = \hat\gamma_k/\hat\gamma_0}, \eqn{\hat\gamma_k=\frac1n\sum_{t}
#' (x_t-\bar x)(x_{t-k}-\bar x)}.
#'
#' @param x serie 
#' @param lag.max retard maximal
#' @return vecteur des autocorrelations (retards 0 a `lag.max`).
#' @export
acf_ts <- function(x, lag.max = 10L) {
  x <- x - mean(x); n <- length(x); c0 <- sum(x^2) / n
  vapply(0:lag.max, function(l) sum(x[(l + 1):n] * x[1:(n - l)]) / n / c0, numeric(1))
}

#' Fonction d'autocorrelation partielle (PACF, recursion de Durbin-Levinson)
#'
#' @param x serie 
#' @param lag.max retard maximal
#' @return vecteur des autocorrelations partielles (retards 1 a `lag.max`).
#' @export
pacf_ts <- function(x, lag.max = 10L) {
  r <- acf_ts(x, lag.max); k <- lag.max; phi <- matrix(0, k, k); phi[1, 1] <- r[2]
  for (m in 2:k) {
    num <- r[m + 1] - sum(phi[m - 1, 1:(m - 1)] * r[m:2])
    den <- 1 - sum(phi[m - 1, 1:(m - 1)] * r[2:m])
    phi[m, m] <- num / den
    for (j in 1:(m - 1)) phi[m, j] <- phi[m - 1, j] - phi[m, m] * phi[m - 1, m - j]
  }
  diag(phi)
}

#' Estimation d'un AR(p) par les equations de Yule-Walker
#'
#' Resout \eqn{R\phi=r} (R matrice de Toeplitz des autocorrelations). Estimateur
#' de moments, coherent et rapide.
#'
#' @param x serie 
#' @param order ordre p
#' @return liste : `ar` (coefficients), `var_pred` (variance d'innovation), `mean`.
#' @export
ar_yw <- function(x, order = 1L) {
  r <- acf_ts(x, order); R <- toeplitz(r[1:order]); phi <- as.numeric(solve(R, r[2:(order + 1)]))
  s2 <- (sum((x - mean(x))^2) / length(x)) * (1 - sum(phi * r[2:(order + 1)]))
  list(ar = phi, var_pred = s2, mean = mean(x))
}

#' Estimation d'un ARMA(p,q) par moindres carres conditionnels (CSS)
#'
#' Minimise \eqn{\sum_t \hat\varepsilon_t^2}, ou \eqn{\hat\varepsilon_t=(x_t-\mu)-
#' \sum \phi_i(x_{t-i}-\mu)-\sum\theta_j\hat\varepsilon_{t-j}} est calcule
#' recursivement (erreurs initiales nulles). Optimisation via `optim`.
#'
#' @param x serie 
#' @param p,q ordres AR et MA
#' @return liste : `ar`, `ma`, `mean`, `sigma2`.
#' @export
arma_css <- function(x, p = 1L, q = 1L) {
  n <- length(x); m <- max(p, q)
  css <- function(par) {
    mu <- par[1]; ar <- par[1 + seq_len(p)]; ma <- par[1 + p + seq_len(q)]
    z <- x - mu; e <- numeric(n)
    for (t in seq_len(n)) {
      ARt <- if (t > p) sum(ar * z[(t - 1):(t - p)]) else 0
      MAt <- if (t > q) sum(ma * e[(t - 1):(t - q)]) else 0
      e[t] <- z[t] - ARt - MAt
    }
    sum(e[(m + 1):n]^2)
  }
  o <- optim(c(mean(x), rep(0.1, p + q)), css, method = "BFGS")
  list(ar = o$par[1 + seq_len(p)], ma = o$par[1 + p + seq_len(q)],
       mean = o$par[1], sigma2 = o$value / (n - m))
}

#' Test portmanteau de Ljung-Box (autocorrelation residuelle)
#'
#' \eqn{Q=n(n+2)\sum_{k=1}^{K}\hat\rho_k^2/(n-k)\sim\chi^2_K} sous absence
#' d'autocorrelation. Utilise sur les residus d'un modele ajuste.
#'
#' @param x serie (ou residus) 
#' @param lag nombre de retards K
#' @param fitdf degres de liberte du modele ajuste (a soustraire).
#' @return liste : `statistic`, `df`, `p_value`.
#' @export
ljung_box <- function(x, lag = 10L, fitdf = 0L) {
  r <- acf_ts(x, lag)[-1]; n <- length(x)
  Q <- n * (n + 2) * sum(r^2 / (n - seq_len(lag)))
  df <- lag - fitdf
  list(statistic = Q, df = df, p_value = pchisq(Q, df, lower.tail = FALSE))
}

#' Test de racine unitaire de Dickey-Fuller augmente (ADF)
#'
#' Regression \eqn{\Delta x_t=\alpha+\beta t+\gamma x_{t-1}+\sum_{i=1}^{k}
#' \delta_i\Delta x_{t-i}+\varepsilon_t} ; la statistique est le t de \eqn{\gamma}.
#' \eqn{H_0:\gamma=0} (racine unitaire, NON stationnaire). Valeurs critiques de
#' Dickey-Fuller (non gaussiennes) : rejet si la stat est TRES negative.
#'
#' @param x serie 
#' @param lags nombre de retards \eqn{k} (defaut : Schwert)
#' @return liste : `statistic`, `lags`.
#' @export
adf_test <- function(x, lags = trunc((length(x) - 1)^(1 / 3))) {
  n <- length(x); dx <- diff(x); k <- lags
  D <- embed(dx, k + 1); dy <- D[, 1]; lagsdx <- D[, -1, drop = FALSE]
  yl <- x[k:(n - 1)][seq_along(dy)]; trend <- ((k + 1):n)[seq_along(dy)]
  X <- cbind(1, trend, yl, lagsdx)
  b <- solve(crossprod(X), crossprod(X, dy)); r <- dy - X %*% b
  s2 <- sum(r^2) / (length(dy) - ncol(X)); V <- s2 * solve(crossprod(X))
  list(statistic = as.numeric(b[3] / sqrt(V[3, 3])), lags = k)
}
