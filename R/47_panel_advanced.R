# =============================================================================
# Module 47 — Panel avance : controle synthetique et panel dynamique
# Implemente les equations de derivations/47_panel_advanced.qmd. R base + quadprog.
# Deux methodes de causalite en panel : le CONTROLE SYNTHETIQUE (construire un
# "clone" pondere de l'unite traitee a partir des controles) et le PANEL
# DYNAMIQUE (estimer un y_{t-1} autoregressif sans le biais de Nickell des EF).
# =============================================================================

#' Controle synthetique (Abadie)
#'
#' Construit un "contrefactuel" de l'unite traitee comme combinaison **convexe**
#' des unites de controle (poids \eqn{w\ge0}, \eqn{\sum w=1}) qui reproduit au
#' mieux la trajectoire PRE-traitement. L'effet est l'ecart post-traitement entre
#' l'unite traitee et son synthetique. Poids par programmation quadratique.
#'
#' @param Y1 trajectoire de l'unite traitee (longueur T).
#' @param Y0 matrice T x J des unites de controle (donor pool).
#' @param pre indices des periodes PRE-traitement.
#' @return liste : `weights`, `synthetic` (trajectoire), `effect` (post-traitement).
#' @export
synthetic_control <- function(Y1, Y0, pre) {
  Y0 <- as.matrix(Y0); J <- ncol(Y0)
  Dmat <- crossprod(Y0[pre, , drop = FALSE]) + diag(1e-8, J)
  dvec <- crossprod(Y0[pre, , drop = FALSE], Y1[pre])
  Amat <- cbind(rep(1, J), diag(J)); bvec <- c(1, rep(0, J))
  w <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
  w <- pmax(w, 0) / sum(pmax(w, 0))
  synth <- as.numeric(Y0 %*% w)
  list(weights = w, synthetic = synth, effect = (Y1 - synth)[-pre])
}

#' Panel dynamique par variables instrumentales (Anderson-Hsiao / Arellano-Bond)
#'
#' Dans \eqn{y_{it}=\rho y_{i,t-1}+\alpha_i+\varepsilon_{it}}, l'estimateur
#' **within** (effets fixes) est BIAISE (\eqn{y_{i,t-1}} demeaне est correle a
#' l'erreur demeanee : **biais de Nickell**). On differencie pour eliminer
#' \eqn{\alpha_i} puis on instrumente \eqn{\Delta y_{i,t-1}} par le **niveau
#' retarde** \eqn{y_{i,t-2}} (valide car non correle a \eqn{\Delta\varepsilon}).
#'
#' @param data data.frame ; @param id,time,y noms des colonnes.
#' @return liste : `rho` (coefficient autoregressif).
#' @export
dynamic_panel_iv <- function(data, id = "id", time = "time", y = "y") {
  d <- data[order(data[[id]], data[[time]]), ]
  ids <- unique(d[[id]]); num <- den <- 0
  for (i in ids) {
    yi <- d[[y]][d[[id]] == i]; Ti <- length(yi)
    if (Ti < 3) next
    for (t in 3:Ti) {
      dy <- yi[t] - yi[t - 1]; dyl <- yi[t - 1] - yi[t - 2]; z <- yi[t - 2]  # instrument = niveau
      num <- num + z * dy; den <- den + z * dyl
    }
  }
  list(rho = num / den)
}

#' Estimateur within (effets fixes) d'un panel dynamique — pour illustrer le biais
#'
#' @param data,id,time,y cf. `dynamic_panel_iv`.
#' @return le coefficient AR estime (biaise : Nickell).
#' @export
dynamic_panel_fe <- function(data, id = "id", time = "time", y = "y") {
  d <- data[order(data[[id]], data[[time]]), ]
  d$ly <- ave(d[[y]], d[[id]], FUN = function(v) c(NA, v[-length(v)]))
  ok <- !is.na(d$ly)
  yv <- d[[y]][ok]; lv <- d$ly[ok]; gid <- d[[id]][ok]
  dem <- function(v, g) v - ave(v, g)
  ytil <- dem(yv, gid); ltil <- dem(lv, gid)
  sum(ltil * ytil) / sum(ltil^2)
}
