# =============================================================================
# BENCHMARK — implémentations « from scratch » vs packages de référence
# -----------------------------------------------------------------------------
# Trois axes : (A) PRÉCISION (les coefficients coïncident-ils ?),
# (B) VITESSE (temps médian), (C) STABILITÉ numérique (mal-conditionnement).
# Objectif : montrer que le code base-R reproduit les références à la tolérance
# machine, à un coût raisonnable, et qu'il est numériquement SAIN (QR vs
# équations normales) — pas seulement « correct sur des cas faciles ».
# =============================================================================

for (f in c("00_linalg", "01_ols", "03_glm_irls", "04_regularisation", "05_iv_2sls"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages({ library(microbenchmark); library(AER); library(ggplot2) })
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) PRÉCISION — écart maximal des coefficients à la référence
# =============================================================================
acc <- list()

## OLS (M1) vs stats::lm
n <- 500; p <- 8
X <- matrix(rnorm(n * p), n, p); y <- as.numeric(X %*% rnorm(p)) + rnorm(n)
d <- data.frame(y = y, X); form <- as.formula(paste("y ~", paste(names(d)[-1], collapse = "+")))
b_ols <- ols_fit(form, d)$coefficients
b_lm  <- coef(lm(form, d))
acc[["OLS (M1) vs lm"]] <- max(abs(b_ols - b_lm))

## GLM logistique (M3, IRLS) vs stats::glm
pr <- 1 / (1 + exp(-(X %*% (rnorm(p) * 0.5)))); yb <- rbinom(n, 1, pr)
db <- data.frame(y = yb, X)
b_irls <- glm_irls(form, db, family = "binomial")$coefficients
b_glm  <- coef(glm(form, db, family = binomial()))
acc[["GLM logit (M3) vs glm"]] <- max(abs(b_irls - b_glm))

## 2SLS (M5) vs AER::ivreg
z <- rnorm(n); xexo <- rnorm(n)
xend <- 0.8 * z + 0.5 * xexo + rnorm(n)                     # endogène (corrélé au bruit)
yy <- 1 + 2 * xend + 1.5 * xexo + rnorm(n)
b_tsls <- tsls_fit(yy, cbind(1, xend, xexo), cbind(1, z, xexo))$coefficients
b_iv   <- coef(ivreg(yy ~ xend + xexo | z + xexo))
acc[["2SLS (M5) vs ivreg"]] <- max(abs(as.numeric(b_tsls) - as.numeric(b_iv)))

cat("=== (A) PRÉCISION : écart max des coefficients à la référence ===\n\n")
for (nm in names(acc))
  cat(sprintf("  %-24s : %.2e  %s\n", nm, acc[[nm]],
              ifelse(acc[[nm]] < 1e-6, "OK (tolérance machine)", "!! écart")))

# =============================================================================
# (B) VITESSE — temps médian (from scratch vs référence)
# =============================================================================
bm <- microbenchmark(
  `OLS M1`   = ols_fit(form, d),
  `OLS lm`   = lm(form, d),
  `GLM M3`   = glm_irls(form, db, family = "binomial"),
  `GLM glm`  = glm(form, db, family = binomial()),
  `2SLS M5`  = tsls_fit(yy, cbind(1, xend, xexo), cbind(1, z, xexo)),
  `2SLS ivreg` = ivreg(yy ~ xend + xexo | z + xexo),
  times = 50L)
sm <- aggregate(time ~ expr, bm, median); sm$ms <- sm$time / 1e6
cat("\n=== (B) VITESSE : temps médian (ms), n =", n, "===\n\n")
for (i in seq_len(nrow(sm))) cat(sprintf("  %-12s %8.3f ms\n", sm$expr[i], sm$ms[i]))

# =============================================================================
# (C) STABILITÉ — matrice de Läuchli : QR (M0) vs équations normales (Cholesky)
# =============================================================================
# Exemple canonique (Läuchli). Avec X = [[1,1,1],[eps,0,0],[0,eps,0],[0,0,eps]],
# X'X a une diagonale 1+eps^2 et un conditionnement ~ 1/eps^2 : former X'X
# (équations normales / Cholesky) CARRE le conditionnement et perd ~la moitié
# des chiffres significatifs ; la factorisation QR travaille sur X directement
# (conditionnement ~ 1/eps) et reste précise. Ici eps=1e-6 : QR ~1/eps=1e6,
# équations normales ~1/eps^2=1e12.
eps <- 1e-6
Xi <- rbind(c(1, 1, 1), diag(eps, 3))                      # 4 x 3
beta_true <- c(1, -2, 3); yi <- as.numeric(Xi %*% beta_true)   # système consistant
b_qr   <- solve_ls_qr(Xi, yi)$coefficients                 # QR (M0)
b_chol <- tryCatch(solve_ls_chol(Xi, yi)$coefficients, error = function(e) rep(NA_real_, 3))
b_lmi  <- tryCatch(coef(lm(yi ~ Xi - 1)), error = function(e) rep(NA_real_, 3))
err <- function(b) if (all(is.finite(b))) max(abs(b - beta_true)) else Inf
cat(sprintf("\n=== (C) STABILITÉ : matrice de Läuchli, conditionnement(X) ~ %.0e ===\n\n", 1 / eps))
cat(sprintf("  QR (M0, solve_ls_qr)      : erreur %.2e\n", err(b_qr)))
cat(sprintf("  Équations normales (Chol) : erreur %.2e\n", err(b_chol)))
cat(sprintf("  lm (référence, QR)        : erreur %.2e\n", err(b_lmi)))
cat("\n=> Les équations normales perdent ~la moitié des chiffres significatifs\n")
cat("   (conditionnement au carré) ; QR — comme lm — reste précis. D'où le choix\n")
cat("   systématique de QR dans tout le projet (Module 0).\n")

# --- Graphique vitesse -------------------------------------------------------
sm$modele <- sub(" .*", "", sm$expr); sm$impl <- ifelse(grepl("M[0-9]", sm$expr), "from scratch", "référence")
gg <- ggplot(sm, aes(modele, ms, fill = impl)) +
  geom_col(position = position_dodge()) +
  labs(title = "Vitesse : from scratch vs référence",
       subtitle = paste0("temps médian (ms), n = ", n),
       x = NULL, y = "temps médian (ms)", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "benchmark_speed.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "benchmark_speed.png"), "\n")
