# =============================================================================
# Module 25 — Différences-de-différences à adoption échelonnée
# Implémente les équations de derivations/25_did.qmd. R base.
# Convention des données : panel ÉQUILIBRÉ, data.frame avec colonnes
#   id (unité), t (période), y (résultat), g (période de PREMIER traitement ;
#   Inf ou 0 pour les jamais-traités). Traitement D_it = 1{ t >= g }.
# =============================================================================

# double démean (effets fixes unité+temps ; forme close, panel équilibré) --------
.tw_demean <- function(v, id, t) v - ave(v, id) - ave(v, t) + mean(v)
.is_never <- function(g) is.infinite(g) | g == 0

#' Estimateur TWFE (two-way fixed effects) de l'effet du traitement
#'
#' Régression de \eqn{y} sur l'indicatrice \eqn{D_{it}} avec effets fixes unité
#' et temps (double démean, FWL — Modules 1, 21). SE groupée par unité.
#'
#' @param data data.frame du panel.
#' @param yname,idname,tname,gname noms des colonnes (résultat, unité, temps,
#'   cohorte de premier traitement ; jamais-traités = Inf ou 0).
#' @return liste : `coef`, `se` (groupée par unité).
#' @export
twfe <- function(data, yname = "y", idname = "id", tname = "t", gname = "g") {
  id <- data[[idname]]; t <- data[[tname]]; g <- data[[gname]]
  D <- as.numeric(!.is_never(g) & t >= g)
  yt <- .tw_demean(data[[yname]], id, t); Dt <- .tw_demean(D, id, t)
  coef <- sum(Dt * yt) / sum(Dt^2)
  r <- yt - coef * Dt
  # SE groupée par unité (sandwich clusterisé, éq. 21.3)
  s2 <- tapply(Dt * r, id, sum); meat <- sum(s2^2)
  se <- sqrt(meat) / sum(Dt^2)
  list(coef = coef, se = se, Dtilde = Dt, D = D)
}

#' Poids de la régression TWFE (de Chaisemartin & D'Haultfœuille 2020)
#'
#' Le coefficient TWFE est \eqn{\sum_{(i,t):D=1} w_{it}\,\tau_{it}} avec
#' \eqn{w_{it}=\tilde D_{it}/\sum \tilde D^2} (\eqn{\tilde D} = traitement
#' résidualisé des effets fixes). Les poids somment à 1 mais certains sont
#' **négatifs** : TWFE n'est pas une moyenne convexe des effets — il peut même
#' être de signe opposé à tous les effets individuels.
#'
#' @param data data.frame du panel ; noms de colonnes comme `twfe`.
#' @inheritParams twfe
#' @return liste : `weights` (sur cellules traitées), `share_negative`, `sum`.
#' @export
twfe_weights <- function(data, yname = "y", idname = "id", tname = "t", gname = "g") {
  tw <- twfe(data, yname, idname, tname, gname)
  trt <- tw$D == 1
  w <- tw$Dtilde[trt] / sum(tw$Dtilde^2)
  list(weights = w, share_negative = mean(w < 0), sum = sum(w))
}

#' ATT groupe-temps (Callaway & Sant'Anna 2021)
#'
#' \eqn{\mathrm{ATT}(g,t)=\bigl[\bar Y_{g,t}-\bar Y_{g,g-1}\bigr]-
#' \bigl[\bar Y_{C,t}-\bar Y_{C,g-1}\bigr]}, base universelle \eqn{g-1}, contrôles
#' = jamais-traités (`"never"`) ou pas-encore-traités (`"notyet"`). Évite les
#' « comparaisons interdites » du TWFE.
#'
#' @param data data.frame du panel ; noms comme `twfe`.
#' @inheritParams twfe
#' @param control "never" (jamais-traités) ou "notyet" (pas encore traités en t).
#' @return data.frame : `g`, `t`, `att`, `n_treat`.
#' @export
att_gt <- function(data, yname = "y", idname = "id", tname = "t", gname = "g",
                   control = c("never", "notyet")) {
  control <- match.arg(control)
  y <- data[[yname]]; id <- data[[idname]]; t <- data[[tname]]; g <- data[[gname]]
  ybar <- function(mask, per) mean(y[mask & t == per])
  gs <- sort(unique(g[!.is_never(g)])); ts <- sort(unique(t))
  out <- list()
  for (gg in gs) for (tt in ts) {
    if (tt == gg - 1) next
    base <- gg - 1
    trg <- g == gg
    ctl <- if (control == "never") .is_never(g) else (.is_never(g) | g > max(tt, base))
    d_tr <- ybar(trg, tt) - ybar(trg, base)
    d_co <- ybar(ctl, tt) - ybar(ctl, base)
    out[[length(out) + 1L]] <- data.frame(g = gg, t = tt, att = d_tr - d_co,
                                           n_treat = sum(trg & t == tt))
  }
  do.call(rbind, out)
}

#' Agrégation des ATT(g,t) (Callaway-Sant'Anna)
#'
#' @param attgt sortie de `att_gt`.
#' @param type "simple" (ATT global, périodes post, pondéré par taille de groupe),
#'   "dynamic" (event-study par ancienneté \eqn{e=t-g}) ou "group" (par cohorte).
#' @return pour "simple" un scalaire ; sinon un data.frame (`e`/`g`, `att`).
#' @export
aggregate_att <- function(attgt, type = c("simple", "dynamic", "group")) {
  type <- match.arg(type)
  post <- attgt[attgt$t >= attgt$g, ]
  if (type == "simple") return(sum(post$att * post$n_treat) / sum(post$n_treat))
  if (type == "group") {
    ag <- tapply(seq_len(nrow(post)), post$g, function(i)
      sum(post$att[i] * post$n_treat[i]) / sum(post$n_treat[i]))
    return(data.frame(g = as.numeric(names(ag)), att = as.numeric(ag)))
  }
  post$e <- post$t - post$g
  ad <- tapply(seq_len(nrow(post)), post$e, function(i)
    sum(post$att[i] * post$n_treat[i]) / sum(post$n_treat[i]))
  data.frame(e = as.numeric(names(ad)), att = as.numeric(ad))
}

#' Event-study de Sun & Abraham (2021) — interaction-weighted
#'
#' Régression saturée \eqn{y} sur les indicatrices cohorte×ancienneté (référence
#' \eqn{e=-1}) avec effets fixes unité+temps, puis agrégation des coefficients par
#' ancienneté \eqn{e} pondérée par les parts de cohorte. Robuste aux effets
#' dynamiques hétérogènes, contrairement au TWFE.
#'
#' @param data data.frame du panel ; noms comme `twfe`.
#' @inheritParams twfe
#' @return liste : `es` (data.frame `e`, `att` — event-study) et `att` (ATT global
#'   post-traitement, pondéré par les tailles de cellule).
#' @export
sunab <- function(data, yname = "y", idname = "id", tname = "t", gname = "g") {
  y <- data[[yname]]; id <- data[[idname]]; t <- data[[tname]]; g <- data[[gname]]
  rel <- ifelse(.is_never(g), NA, t - g)
  cells <- unique(data.frame(g = g, rel = rel)[!is.na(rel) & rel != -1, ])
  cells <- cells[order(cells$g, cells$rel), ]
  X <- vapply(seq_len(nrow(cells)), function(k)
    as.numeric(!is.na(rel) & g == cells$g[k] & rel == cells$rel[k]), numeric(length(y)))
  yt <- .tw_demean(y, id, t); Xt <- apply(X, 2, .tw_demean, id = id, t = t)
  bb <- solve(crossprod(Xt), crossprod(Xt, yt))[, 1]
  ncell <- vapply(seq_len(nrow(cells)), function(k)
    sum(!is.na(rel) & g == cells$g[k] & rel == cells$rel[k]), numeric(1))
  evs <- sort(unique(cells$rel))
  es <- vapply(evs, function(e) { k <- cells$rel == e
    sum(ncell[k] * bb[k]) / sum(ncell[k]) }, numeric(1))
  post <- cells$rel >= 0
  list(es = data.frame(e = evs, att = es),
       att = sum(ncell[post] * bb[post]) / sum(ncell[post]))
}
