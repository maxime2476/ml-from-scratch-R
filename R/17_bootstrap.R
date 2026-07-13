# =============================================================================
# Module 17 — Bootstrap
# Implémente les équations de derivations/17_bootstrap.qmd. R base + Module 0.
# =============================================================================

#' Bootstrap non paramétrique d'une statistique
#'
#' Rééchantillonne les données avec remise et applique `stat_fn` (principe
#' plug-in). Renvoie les répliques et l'erreur standard bootstrap (éq. 17.1).
#'
#' @param data vecteur, matrice ou data.frame (les lignes sont rééchantillonnées).
#' @param stat_fn fonction `data -> statistique scalaire`.
#' @param R nombre de rééchantillons.
#' @param seed graine.
#' @return objet `bootstrap` : `t0`, `replicates`, `se`, `data`, `stat_fn`, `n`.
#' @export
bootstrap <- function(data, stat_fn, R = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  is_df <- is.data.frame(data) || is.matrix(data)
  n <- if (is_df) nrow(data) else length(data)
  reps <- numeric(R)
  for (b in seq_len(R)) {
    idx <- sample.int(n, n, replace = TRUE)
    reps[b] <- stat_fn(if (is_df) data[idx, , drop = FALSE] else data[idx])
  }
  structure(list(t0 = stat_fn(data), replicates = reps, se = sd(reps),
                 data = data, stat_fn = stat_fn, n = n, is_df = is_df),
            class = "bootstrap")
}

# Accélération BCa par jackknife (éq. 17.2).
.bca_ci <- function(bt, level) {
  reps <- bt$replicates; t0 <- bt$t0
  z0 <- qnorm(mean(reps < t0) + 0.5 * mean(reps == t0))
  jack <- vapply(seq_len(bt$n), function(i)
    bt$stat_fn(if (bt$is_df) bt$data[-i, , drop = FALSE] else bt$data[-i]), numeric(1))
  u <- mean(jack) - jack
  a <- sum(u^3) / (6 * (sum(u^2))^1.5)
  alpha <- (1 - level) / 2; za <- qnorm(c(alpha, 1 - alpha))
  adj <- pnorm(z0 + (z0 + za) / (1 - a * (z0 + za)))
  unname(quantile(reps, adj, names = FALSE, type = 6))
}

#' Intervalles de confiance bootstrap
#'
#' Types : "percentile", "basic" (pivotal), "normal", "bca" (bias-corrected
#' accelerated, éq. 17.2).
#'
#' @param bt objet `bootstrap`.
#' @param level niveau de confiance.
#' @param type "percentile", "basic", "normal" ou "bca".
#' @return vecteur (borne inf, borne sup).
#' @export
boot_ci <- function(bt, level = 0.95, type = c("percentile", "basic", "normal", "bca")) {
  type <- match.arg(type); a <- 1 - level; reps <- bt$replicates; t0 <- bt$t0
  qs <- c(a / 2, 1 - a / 2)
  switch(type,
    percentile = unname(quantile(reps, qs, names = FALSE, type = 6)),
    basic = { q <- unname(quantile(reps, rev(qs), names = FALSE, type = 6)); 2 * t0 - q },
    normal = t0 + qnorm(qs) * bt$se,
    bca = .bca_ci(bt, level))
}

#' Bootstrap d'une régression linéaire (pairs ou résidus)
#'
#' "pairs" rééchantillonne les couples \eqn{(x_i,y_i)} (robuste à
#' l'hétéroscédasticité ; analogue du sandwich) ; "residual" rééchantillonne les
#' résidus (suppose des erreurs i.i.d.).
#'
#' @param formula formule.
#' @param data data.frame.
#' @param R nombre de rééchantillons.
#' @param method "pairs" ou "residual".
#' @param seed graine.
#' @return objet `boot_lm` : `t0` (coefficients), `replicates` (R x p), `se`.
#' @export
boot_lm <- function(formula, data, R = 2000L, method = c("pairs", "residual"),
                    seed = NULL) {
  method <- match.arg(method); if (!is.null(seed)) set.seed(seed)
  mf <- model.frame(formula, data)
  y <- as.numeric(model.response(mf)); X <- model.matrix(attr(mf, "terms"), mf)
  n <- nrow(X); p <- ncol(X)
  b0 <- solve_ls_qr(X, y)$coefficients; names(b0) <- colnames(X)
  fitted <- as.numeric(X %*% b0); resid <- y - fitted
  reps <- matrix(NA_real_, R, p, dimnames = list(NULL, colnames(X)))
  for (r in seq_len(R)) {
    if (method == "pairs") {
      idx <- sample.int(n, n, replace = TRUE)
      reps[r, ] <- solve_ls_qr(X[idx, , drop = FALSE], y[idx])$coefficients
    } else {
      ystar <- fitted + sample(resid, n, replace = TRUE)
      reps[r, ] <- solve_ls_qr(X, ystar)$coefficients
    }
  }
  structure(list(t0 = b0, replicates = reps, se = apply(reps, 2, sd),
                 method = method), class = "boot_lm")
}
