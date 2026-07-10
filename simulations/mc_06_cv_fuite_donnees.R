# =============================================================================
# Monte Carlo — Module 6 : fuite de données en validation croisée
# Démontre le BIAIS OPTIMISTE induit quand un prétraitement voit les données de
# test. Deux niveaux :
#  (A) Standardisation AVANT vs APRÈS le découpage (fuite douce : ne touche que
#      les moments de X, biais faible).
#  (B) Sélection de variables supervisée AVANT vs DANS les plis (fuite via les
#      LABELS : biais énorme — l'exemple classique ESL §7.10.2). Avec y
#      INDÉPENDANT de X, la « mauvaise » CV fabrique une compétence prédictive
#      fictive ; la « bonne » CV révèle l'absence de signal.
# =============================================================================

for (f in c("00_linalg", "06_validation"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_folds <- function(n, K) sample(rep_len(seq_len(K), n))

# =============================================================================
# (A) Standardisation avant vs après le split
# =============================================================================
# DGP avec un vrai signal ; on compare l'estimation CV de l'erreur de test à
# l'erreur de test VRAIE (grand échantillon indépendant).
nA <- 80; pA <- 20; K <- 10; R_A <- 400
biasA <- data.frame(fuite = numeric(R_A), correct = numeric(R_A), vrai = numeric(R_A))

for (r in seq_len(R_A)) {
  X <- matrix(rnorm(nA * pA), nA, pA)
  beta <- c(rep(1, 5), rep(0, pA - 5))
  y <- as.numeric(X %*% beta + rnorm(nA, sd = 3))
  folds <- make_folds(nA, K)

  # -- CV avec FUITE : on standardise X et y sur TOUT l'échantillon, puis CV --
  mx <- colMeans(X); sx <- apply(X, 2, sd)
  Xf <- scale(X, mx, sx); yf <- y - mean(y)
  errf <- numeric(nA)
  for (k in seq_len(K)) {
    te <- which(folds == k); tr <- which(folds != k)
    b <- solve_ls_qr(cbind(1, Xf[tr, ]), yf[tr])$coefficients
    errf[te] <- (yf[te] - cbind(1, Xf[te, ]) %*% b)^2
  }

  # -- CV CORRECTE : standardisation calculée sur le TRAIN de chaque pli --
  errc <- numeric(nA)
  for (k in seq_len(K)) {
    te <- which(folds == k); tr <- which(folds != k)
    mxt <- colMeans(X[tr, ]); sxt <- apply(X[tr, ], 2, sd); myt <- mean(y[tr])
    Xtr <- scale(X[tr, ], mxt, sxt); Xte <- scale(X[te, ], mxt, sxt)
    b <- solve_ls_qr(cbind(1, Xtr), y[tr] - myt)$coefficients
    errc[te] <- ((y[te] - myt) - cbind(1, Xte) %*% b)^2
  }

  # -- Erreur de test VRAIE (grand échantillon indépendant) --
  Xtest <- matrix(rnorm(5000 * pA), 5000, pA)
  ytest <- as.numeric(Xtest %*% beta + rnorm(5000, sd = 3))
  b_full <- solve_ls_qr(cbind(1, scale(X, mx, sx)), y - mean(y))$coefficients
  pred <- mean(y) + cbind(1, scale(Xtest, mx, sx)) %*% b_full
  biasA[r, ] <- c(mean(errf), mean(errc), mean((ytest - pred)^2))
}

cat("=== (A) Standardisation avant/après le split (erreur quadratique) ===\n")
cat(sprintf("CV avec fuite (standardisation sur tout)   : %.3f\n", mean(biasA$fuite)))
cat(sprintf("CV correcte  (standardisation sur le train): %.3f\n", mean(biasA$correct)))
cat(sprintf("Erreur de test vraie                        : %.3f\n", mean(biasA$vrai)))
cat(sprintf("EFFET DE LA FUITE (fuite - correcte)        : %+.4f\n",
            mean(biasA$fuite - biasA$correct)))
cat("=> quasi nul : la standardisation NON supervisée de X ne fuit presque pas.\n")

# =============================================================================
# (B) Sélection de variables supervisée : la fuite qui compte
# =============================================================================
# y INDÉPENDANT de X (aucun signal) ; on sélectionne les k variables les plus
# corrélées à y, puis on ajuste. Erreur de test vraie = Var(y).
nB <- 50; pB <- 1000; ksel <- 10; KB <- 5; R_B <- 300
resB <- data.frame(fuite = numeric(R_B), correct = numeric(R_B))

for (r in seq_len(R_B)) {
  X <- matrix(rnorm(nB * pB), nB, pB)
  y <- rnorm(nB)                                   # AUCUN lien avec X
  folds <- make_folds(nB, KB)

  # -- MAUVAISE façon : sélectionner sur TOUT l'échantillon, puis CV --
  cors <- abs(as.numeric(cor(X, y)))
  sel <- order(cors, decreasing = TRUE)[seq_len(ksel)]
  errf <- numeric(nB)
  for (k in seq_len(KB)) {
    te <- which(folds == k); tr <- which(folds != k)
    b <- solve_ls_qr(cbind(1, X[tr, sel]), y[tr])$coefficients
    errf[te] <- (y[te] - cbind(1, X[te, sel]) %*% b)^2
  }

  # -- BONNE façon : sélectionner DANS chaque pli (sur le train seulement) --
  errc <- numeric(nB)
  for (k in seq_len(KB)) {
    te <- which(folds == k); tr <- which(folds != k)
    cors_tr <- abs(as.numeric(cor(X[tr, ], y[tr])))
    sel_tr <- order(cors_tr, decreasing = TRUE)[seq_len(ksel)]
    b <- solve_ls_qr(cbind(1, X[tr, sel_tr]), y[tr])$coefficients
    errc[te] <- (y[te] - cbind(1, X[te, sel_tr]) %*% b)^2
  }
  resB[r, ] <- c(mean(errf), mean(errc))
}

vary <- 1  # Var(y) = 1 (erreur de test vraie, aucun signal)
cat("\n=== (B) Sélection supervisée avant vs dans les plis (y sans lien avec X) ===\n")
cat(sprintf("Var(y) (erreur de test vraie) : %.3f\n", vary))
cat(sprintf("CV avec fuite (sélection avant split)  : %.3f  -> pseudo-R² = %+.2f\n",
            mean(resB$fuite), 1 - mean(resB$fuite) / vary))
cat(sprintf("CV correcte (sélection dans les plis)  : %.3f  -> pseudo-R² = %+.2f\n",
            mean(resB$correct), 1 - mean(resB$correct) / vary))
cat("La fuite par sélection supervisée fabrique une compétence prédictive fictive.\n")

# ---- Graphique : distributions des erreurs CV (partie B) --------------------
dfB <- rbind(data.frame(err = resB$fuite, methode = "Fuite : sélection avant le split"),
             data.frame(err = resB$correct, methode = "Correct : sélection dans les plis"))
ggB <- ggplot(dfB, aes(err, fill = methode)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = vary, linetype = "dashed") +
  annotate("text", x = vary, y = 0, label = "Var(y) = erreur vraie", vjust = -0.5, hjust = 1.05, size = 3) +
  labs(title = "Fuite de données : sélection de variables supervisée (y sans signal)",
       subtitle = "La CV avec fuite sous-estime largement l'erreur (fausse compétence) ; la CV correcte retrouve Var(y)",
       x = "erreur CV estimée", y = "densité", fill = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_06_cv_fuite.png"), ggB, width = 8, height = 5, dpi = 120)

cat("\nGraphique -> ", file.path(out_dir, "mc_06_cv_fuite.png"), "\n")
