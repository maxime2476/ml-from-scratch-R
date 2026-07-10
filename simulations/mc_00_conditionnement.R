# =============================================================================
# Monte Carlo — Module 0 : conditionnement et perte de précision
# Illustre numériquement la Proposition 0.1 (kappa(X'X) = kappa(X)^2) : sur un
# design de conditionnement kappa croissant contrôlé, la résolution par
# équations normales (Cholesky / inverse) subit kappa^2 tandis que QR et SVD
# préservent kappa. On mesure l'erreur relative de reconstruction de beta_vrai.
#
# DGP : X = U diag(sigma) V^T avec sigma de 1 à 1/kappa (SVD imposée) ; beta
# vrai fixé ; y = X beta (SANS bruit, la solution exacte existe et vaut beta).
# Toute erreur observée est donc purement numérique.
# =============================================================================

source("R/00_linalg.R")   # exécuter depuis la racine du projet (ml-from-scratch-R/)
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Générateur de design à conditionnement imposé --------------------------
make_illcond_X <- function(n, p, kappa) {
  U <- qr.Q(qr(matrix(rnorm(n * p), n, p)))   # n x p, colonnes orthonormées
  V <- qr.Q(qr(matrix(rnorm(p * p), p, p)))   # p x p orthogonale
  sigma <- exp(seq(0, -log(kappa), length.out = p))  # sigma_1=1 ... sigma_p=1/kappa
  (U * rep(sigma, each = n)) %*% t(V)         # U diag(sigma) V^T
}

# ---- Résolution par équations normales explicites (base R) ------------------
# Renvoie NULL si solve() déclare le système singulier (= divergence numérique).
solve_ls_normal <- function(X, y)
  tryCatch(as.numeric(solve(crossprod(X), crossprod(X, y))), error = function(e) NULL)

# ---- Balayage du conditionnement --------------------------------------------
n <- 200; p <- 8
kappas <- 10^seq(1, 9, by = 0.5)
n_rep  <- 300               # illustration numérique (pas d'aléa inférentiel ici)
beta_true <- rep(1, p)

grid <- expand.grid(kappa = kappas, methode = c("QR", "Cholesky", "Normale (solve)", "SVD"),
                    stringsAsFactors = FALSE)
grid$err <- NA_real_

err_rel <- function(bhat) sqrt(sum((bhat - beta_true)^2)) / sqrt(sum(beta_true^2))

for (kap in kappas) {
  acc <- c(QR = 0, Cholesky = 0, `Normale (solve)` = 0, SVD = 0)
  n_ok <- c(Cholesky = 0L, `Normale (solve)` = 0L)
  for (r in seq_len(n_rep)) {
    X <- make_illcond_X(n, p, kap)
    y <- as.numeric(X %*% beta_true)
    acc["QR"]  <- acc["QR"]  + err_rel(solve_ls_qr(X, y)$coefficients)
    acc["SVD"] <- acc["SVD"] + err_rel(solve_ls_svd(X, y)$coefficients)
    # Équations normales : Cholesky et solve() peuvent échouer (X'X numériquement
    # non SPD) dès que kappa^2 dépasse la précision -> c'est la divergence.
    bn <- solve_ls_normal(X, y)
    if (!is.null(bn)) { acc["Normale (solve)"] <- acc["Normale (solve)"] + err_rel(bn)
                        n_ok["Normale (solve)"] <- n_ok["Normale (solve)"] + 1L }
    bc <- tryCatch(solve_ls_chol(X, y)$coefficients, error = function(e) NULL)
    if (!is.null(bc)) { acc["Cholesky"] <- acc["Cholesky"] + err_rel(bc)
                        n_ok["Cholesky"] <- n_ok["Cholesky"] + 1L }
  }
  grid$err[grid$kappa == kap & grid$methode == "QR"]  <- acc["QR"]  / n_rep
  grid$err[grid$kappa == kap & grid$methode == "SVD"] <- acc["SVD"] / n_rep
  grid$err[grid$kappa == kap & grid$methode == "Cholesky"] <-
    if (n_ok["Cholesky"] > 0) acc["Cholesky"] / n_ok["Cholesky"] else NA_real_
  grid$err[grid$kappa == kap & grid$methode == "Normale (solve)"] <-
    if (n_ok["Normale (solve)"] > 0) acc["Normale (solve)"] / n_ok["Normale (solve)"] else NA_real_
}

# ---- Tableau récapitulatif ---------------------------------------------------
tab <- reshape(grid, idvar = "kappa", timevar = "methode", direction = "wide")
names(tab) <- sub("err\\.", "", names(tab))
cat("\n=== Erreur relative moyenne de reconstruction de beta (", n_rep, "rép.) ===\n")
print(format(tab, digits = 3, scientific = TRUE), row.names = FALSE)

eps <- .Machine$double.eps
cat("\nPrécision machine eps =", eps, "\n")
cat("Prédiction théorique : erreur QR/SVD ~ kappa*eps ; erreur Normale/Chol ~ kappa^2*eps.\n")

# ---- Graphique ---------------------------------------------------------------
grid$type <- ifelse(grid$methode %in% c("QR", "SVD"), "Sur X (kappa)", "Sur X'X (kappa^2)")
ref <- data.frame(
  kappa = rep(kappas, 2),
  err   = c(kappas * eps, kappas^2 * eps),
  ref   = rep(c("kappa * eps", "kappa^2 * eps"), each = length(kappas))
)

p_plot <- ggplot(grid, aes(kappa, pmax(err, 1e-18), colour = methode, linetype = type)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  geom_line(data = ref, aes(kappa, err, group = ref), colour = "grey40",
            linetype = "dotted", inherit.aes = FALSE) +
  annotate("text", x = 1e5, y = 1e5 * eps * 3, label = "kappa * eps",
           colour = "grey30", size = 3) +
  annotate("text", x = 3e3, y = (3e3)^2 * eps * 3, label = "kappa^2 * eps",
           colour = "grey30", size = 3) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Perte de précision selon le solveur",
       subtitle = expression("QR/SVD suivent "*kappa*" ; équations normales suivent "*kappa^2),
       x = expression(kappa[2](X)), y = "Erreur relative sur "~hat(beta),
       colour = "Méthode", linetype = "Sensibilité") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "mc_00_conditionnement.png"), p_plot,
       width = 8, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "mc_00_conditionnement.png"), "\n")

# ---- Interprétation ----------------------------------------------------------
# - QR et SVD travaillent sur X : erreur ~ kappa*eps, ils restent utilisables
#   jusqu'à kappa ~ 1e15.
# - Les équations normales (Cholesky, ou solve(X'X)) forment X'X et subissent
#   kappa^2 : dès kappa ~ 1e8, l'erreur atteint l'ordre 1 (aucun chiffre fiable),
#   et Cholesky finit par échouer (X'X numériquement non définie positive).
# - Conclusion pratique : préférer QR (forme fermée MCO) ; réserver les
#   équations normales aux problèmes bien conditionnés (ridge, cf. Module 4).
