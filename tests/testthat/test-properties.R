# =============================================================================
# Tests basés sur les PROPRIÉTÉS (property-based / métamorphiques)
# -----------------------------------------------------------------------------
# Au-delà de la validation ponctuelle : des INVARIANTS mathématiques qui doivent
# tenir pour TOUTE entrée. On les éprouve par « fuzzing » léger (plusieurs jeux
# aléatoires). Ces tests attrapent des bugs que les cas fixes manquent.
# =============================================================================

set.seed(20260714)
gen_reg <- function(n = 60, p = 3) { X <- matrix(rnorm(n * p), n, p)
  list(X = X, y = as.numeric(X %*% rnorm(p)) + rnorm(n)) }

test_that("OLS : résidus orthogonaux aux régresseurs (X'e = 0) pour toute entrée", {
  for (s in 1:8) { set.seed(s); d <- gen_reg()
    fit <- ols_fit(y ~ ., data.frame(y = d$y, d$X))
    Xd <- model.matrix(y ~ ., data.frame(y = d$y, d$X))
    expect_lt(max(abs(crossprod(Xd, fit$residuals))), 1e-8)
  }
})

test_that("OLS : équivariance d'échelle (rescaler x_j -> coef /c, ajustement inchangé)", {
  for (s in 1:8) { set.seed(s); d <- gen_reg(); c <- runif(1, 0.2, 5)
    f1 <- ols_fit(y ~ ., data.frame(y = d$y, d$X))
    X2 <- d$X; X2[, 1] <- X2[, 1] * c
    f2 <- ols_fit(y ~ ., data.frame(y = d$y, X2))
    expect_equal(unname(f2$coefficients[2]), unname(f1$coefficients[2] / c), tolerance = 1e-8)
    expect_equal(f2$fitted, f1$fitted, tolerance = 1e-8)   # invariance de l'ajustement
  }
})

test_that("Fonction d'influence : E[IC] = 0 exactement (M24) pour toute entrée", {
  for (s in 1:8) { set.seed(s); d <- gen_reg()
    X <- cbind(1, d$X)
    expect_lt(max(abs(colMeans(influence_ols(X, d$y)$ic))), 1e-8)
  }
})

test_that("Ridge : rétrécissement monotone (||beta(lambda)|| décroît avec lambda)", {
  for (s in 1:6) { set.seed(s); d <- gen_reg()
    lams <- c(0.01, 1, 10, 100)
    norms <- sapply(lams, function(l) sqrt(sum(ridge_fit(d$X, d$y, lambda = l, intercept = TRUE)$beta^2)))
    expect_true(all(diff(norms) <= 1e-8))               # décroissant
  }
})

test_that("Lasso : la parcimonie croît avec lambda (nb de zéros non décroissant)", {
  for (s in 1:6) { set.seed(s); d <- gen_reg(p = 8)
    lams <- c(0.05, 0.2, 0.5, 1) * max(abs(crossprod(scale(d$X), d$y - mean(d$y))))
    nz <- sapply(lams, function(l) sum(lasso_fit(d$X, d$y, lambda = l)$beta != 0))
    expect_true(all(diff(nz) <= 0))                     # nb de non-nuls non croissant
  }
})

test_that("GLM logistique : probabilités ajustées dans (0,1)", {
  for (s in 1:6) { set.seed(s); n <- 80; x <- rnorm(n); y <- rbinom(n, 1, plogis(x))
    fit <- glm_irls(y ~ x, data.frame(y = y, x = x), family = "binomial")
    expect_true(all(fit$fitted > 0 & fit$fitted < 1))
  }
})

test_that("Noyau RBF : matrice de Gram semi-définie positive (M27)", {
  for (s in 1:6) { set.seed(s); X <- matrix(rnorm(40), 20, 2)
    K <- rbf_kernel(X, X, lengthscale = runif(1, 0.5, 2), variance = runif(1, 0.5, 2))
    expect_gt(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  }
})

test_that("Autodiff : linéarité du gradient grad(a f + b g) = a grad f + b g (M28)", {
  for (s in 1:6) { set.seed(s); x <- rnorm(3); a <- rnorm(1); b <- rnorm(1)
    gf <- ad_grad(function(v) sum(exp(v)), x)
    gg <- ad_grad(function(v) sum(v * v), x)
    gh <- ad_grad(function(v) a * sum(exp(v)) + b * sum(v * v), x)
    expect_equal(gh, a * gf + b * gg, tolerance = 1e-8)
  }
})

test_that("DiD (métamorphique) : sous effet CONSTANT homogène, TWFE/CS/SA non biaisés (M25)", {
  # Invariant : sans dynamique ni hétérogénéité, aucune « comparaison interdite »
  # n'est nuisible -> les TROIS estimateurs sont non biaisés pour tau.
  set.seed(1); N <- 150; Tt <- 6; tau <- 2; R <- 40
  b_tw <- b_cs <- b_sa <- numeric(R)
  for (r in seq_len(R)) {
    id <- rep(1:N, each = Tt); t <- rep(1:Tt, N)
    g <- rep(sample(c(3, 5, Inf), N, replace = TRUE), each = Tt)
    ai <- rep(rnorm(N), each = Tt); gt <- rep(0.2 * (1:Tt), N)
    y <- ai + gt + tau * as.integer(t >= g) + rnorm(length(id))
    d <- data.frame(id, t, g, y)
    b_tw[r] <- twfe(d)$coef
    b_cs[r] <- aggregate_att(att_gt(d, control = "never"), "simple")
    b_sa[r] <- sunab(d)$att
  }
  expect_lt(abs(mean(b_tw) - tau), 0.06)                # biais moyen ~ 0
  expect_lt(abs(mean(b_cs) - tau), 0.06)
  expect_lt(abs(mean(b_sa) - tau), 0.06)
})

test_that("Jackknife et bootstrap : variances positives (M17, M24)", {
  for (s in 1:6) { set.seed(s); y <- rnorm(50)
    expect_gte(jackknife(data.frame(y = y), function(z) mean(z$y))$var, 0)
    expect_gte(var(bootstrap(data.frame(y = y), function(z) mean(z$y), R = 200)$replicates), 0)
  }
})
