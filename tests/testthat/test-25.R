# Tests — Module 25 (DiD échelonné). Références : fixest (TWFE, Sun-Abraham),
# did (Callaway-Sant'Anna). Panel équilibré, adoption échelonnée, effets dynamiques.

make_did <- function(N = 300, Tt = 8, seed = 1) {
  set.seed(seed); cohorts <- c(3, 5, 7, Inf)
  id <- rep(1:N, each = Tt); t <- rep(1:Tt, N)
  g <- rep(sample(cohorts, N, replace = TRUE), each = Tt)
  ai <- rep(rnorm(N), each = Tt); gt <- rep(0.3 * (1:Tt), N)
  D <- as.integer(t >= g); expo <- pmax(t - g + 1, 0); expo[is.infinite(g)] <- 0
  tau <- ifelse(D == 1, 1 + 0.5 * expo + ifelse(g == 3, 1, 0), 0)   # dynamique + hétéro
  y <- ai + gt + tau + rnorm(length(id))
  list(data = data.frame(id, t, g, y), tau = tau, D = D)
}

test_that("TWFE = fixest::feols à 1e-8", {
  skip_if_not_installed("fixest")
  m <- make_did(); d <- m$data
  tw <- twfe(d)
  df <- transform(d, Dd = as.integer(!is.infinite(g) & t >= g))
  bf <- as.numeric(coef(fixest::feols(y ~ Dd | id + t, df, notes = FALSE))["Dd"])
  expect_lt(abs(tw$coef - bf), 1e-8)
})

test_that("poids TWFE : somment à 1, une part est négative (dCDH)", {
  w <- twfe_weights(make_did()$data)
  expect_equal(w$sum, 1, tolerance = 1e-8)
  expect_gt(w$share_negative, 0)                # comparaisons interdites
})

test_that("ATT(g,t) = did::att_gt à 1e-8", {
  skip_if_not_installed("did")
  m <- make_did(); d <- m$data
  ah <- att_gt(d, control = "never"); ah <- ah[ah$t >= ah$g, ]
  dcs <- transform(d, gcs = ifelse(is.infinite(g), 0, g))
  cs <- did::att_gt(yname = "y", tname = "t", idname = "id", gname = "gcs", data = dcs,
                    control_group = "nevertreated", base_period = "universal",
                    est_method = "reg", bstrap = FALSE)
  csd <- data.frame(g = cs$group, t = cs$t, att = cs$att)
  mg <- merge(ah, csd[csd$t >= csd$g, ], by = c("g", "t"))
  expect_lt(max(abs(mg$att.x - mg$att.y)), 1e-8)
})

test_that("Sun-Abraham = fixest::sunab à 1e-8", {
  skip_if_not_installed("fixest")
  m <- make_did(); d <- m$data
  sa <- sunab(d)
  df <- transform(d, gfix = ifelse(is.infinite(g), 10000, g))
  suppressMessages(require(fixest))            # sunab() est un token spécial de la formule fixest
  saf <- summary(fixest::feols(y ~ sunab(gfix, t) | id + t, df, notes = FALSE),
                 agg = "att")$coeftable[1, 1]
  expect_lt(abs(sa$att - saf), 1e-8)
})

test_that("TWFE est biaisé, Callaway-Sant'Anna et Sun-Abraham récupèrent l'ATT vrai", {
  m <- make_did(N = 800); d <- m$data
  att_true <- mean(m$tau[m$D == 1])            # ATT vrai (moyenne des effets traités)
  tw <- twfe(d)$coef
  cs <- aggregate_att(att_gt(d, control = "never"), "simple")
  sa <- sunab(d)$att
  expect_gt(abs(tw - att_true), 0.3)           # TWFE nettement biaisé
  expect_lt(abs(cs - att_true), 0.15)          # CS ~ vrai
  expect_lt(abs(sa - att_true), 0.15)          # SA ~ vrai
})
