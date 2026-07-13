# =============================================================================
# Outils de rigueur pour les études de simulation Monte Carlo
# (Morris, White & Crowther 2019, "Using simulation studies to evaluate
# statistical methods"). Toute quantité estimée par simulation a elle-même une
# ERREUR MONTE CARLO : on la rapporte systématiquement.
# =============================================================================

#' Erreur Monte Carlo de la moyenne d'échantillon
#'
#' \eqn{\mathrm{MCSE}(\bar x) = \hat\sigma/\sqrt R}. C'est l'incertitude due au
#' nombre fini R de réplications (à ne pas confondre avec l'écart-type de
#' l'estimateur étudié).
#'
#' @param x vecteur des R valeurs simulées.
#' @return l'erreur Monte Carlo de leur moyenne.
#' @export
mc_se <- function(x) sd(x) / sqrt(length(x))

#' Résumé d'une étude de simulation : biais, RMSE, variance, avec erreurs MC
#'
#' @param estimates vecteur des R estimations \eqn{\hat\theta}.
#' @param truth valeur vraie \eqn{\theta}.
#' @return liste : `R`, `mean`, `bias` (+ `bias_se`), `rmse` (+ `rmse_se`),
#'   `variance`, `empirical_se`.
#' @export
mc_summary <- function(estimates, truth) {
  R <- length(estimates)
  bias <- mean(estimates) - truth
  sq <- (estimates - truth)^2
  mse <- mean(sq); rmse <- sqrt(mse)
  list(R = R, mean = mean(estimates),
       bias = bias, bias_se = mc_se(estimates),                 # MCSE du biais
       rmse = rmse, rmse_se = mc_se(sq) / (2 * rmse),           # delta-méthode
       variance = var(estimates), empirical_se = sd(estimates))
}

#' Couverture empirique d'IC avec son erreur Monte Carlo (binomiale)
#'
#' Une couverture estimée sur R réplications est une proportion : son erreur MC
#' est \eqn{\sqrt{\hat p(1-\hat p)/R}}. On peut ainsi juger si l'écart à la valeur
#' nominale est significatif.
#'
#' @param covered vecteur logique (l'IC contenait la vraie valeur).
#' @param nominal niveau nominal (défaut 0.95).
#' @return liste : `coverage`, `se`, `ci` (IC de la couverture), `R`, `nominal`,
#'   `nominal_ok` (TRUE si l'IC contient le niveau nominal).
#' @export
coverage_mc <- function(covered, nominal = 0.95) {
  R <- length(covered); p <- mean(covered)
  se <- sqrt(p * (1 - p) / R)
  ci <- p + c(-1, 1) * qnorm(0.975) * se
  list(coverage = p, se = se, ci = ci, R = R, nominal = nominal,
       nominal_ok = ci[1] <= nominal && nominal <= ci[2])
}

#' Taux de rejet (taille/puissance) avec erreur Monte Carlo
#'
#' @param rejected vecteur logique des rejets.
#' @param nominal niveau nominal du test (défaut 0.05, pour la taille).
#' @return liste : `rate`, `se`, `ci`, `R`, `nominal_ok`.
#' @export
reject_mc <- function(rejected, nominal = 0.05) {
  R <- length(rejected); p <- mean(rejected)
  se <- sqrt(p * (1 - p) / R)
  ci <- p + c(-1, 1) * qnorm(0.975) * se
  list(rate = p, se = se, ci = ci, R = R, nominal = nominal,
       nominal_ok = ci[1] <= nominal && nominal <= ci[2])
}

#' Étude de convergence : biais/RMSE et diagnostics de taux selon n
#'
#' Pour chaque taille n, exécute R réplications de `sim_fn(n)` (qui renvoie un
#' \eqn{\hat\theta}) et calcule biais, RMSE (avec erreurs MC), ainsi que les
#' diagnostics de **taux** \eqn{\sqrt n\,\text{biais}} (doit rester borné si le
#' biais est \eqn{O(1/\sqrt n)} ou mieux) et \eqn{\sqrt n\,\hat{sd}} (doit se
#' stabiliser vers l'écart-type asymptotique si l'estimateur est \eqn{\sqrt n}-
#' consistant).
#'
#' @param sim_fn fonction `n -> theta_hat` (une réplication).
#' @param ns vecteur des tailles d'échantillon.
#' @param R réplications par taille.
#' @param truth valeur vraie.
#' @param seed graine.
#' @return data.frame : `n`, `bias`, `bias_se`, `rmse`, `rmse_se`,
#'   `sqrtn_bias`, `sqrtn_sd`.
#' @export
convergence_study <- function(sim_fn, ns, R, truth, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  do.call(rbind, lapply(ns, function(n) {
    est <- vapply(seq_len(R), function(r) sim_fn(n), numeric(1))
    s <- mc_summary(est, truth)
    data.frame(n = n, bias = s$bias, bias_se = s$bias_se,
               rmse = s$rmse, rmse_se = s$rmse_se,
               sqrtn_bias = sqrt(n) * s$bias,
               sqrtn_sd = sqrt(n) * s$empirical_se)
  }))
}

#' Pente log-log du RMSE en fonction de n (taux de convergence empirique)
#'
#' Un estimateur \eqn{\sqrt n}-consistant a un RMSE \eqn{\propto n^{-1/2}}, donc
#' une pente \eqn{\approx -0.5} en échelle log-log.
#'
#' @param conv data.frame issu de `convergence_study`.
#' @return la pente estimée (idéalement ~ -0.5).
#' @export
rmse_rate <- function(conv) unname(coef(lm(log(rmse) ~ log(n), data = conv))[2])
