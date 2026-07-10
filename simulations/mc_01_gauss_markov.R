# =============================================================================
# Monte Carlo — Module 1 : Gauss-Markov, couverture des IC, loi de t
# Vérifie sur un DGP à paramètres vrais connus :
#  (A) non-biais de beta_hat (moyenne des estimations ~ beta_vrai) ;
#  (B) BLUE : un autre estimateur linéaire sans biais a une variance >= OLS ;
#  (C) couverture empirique des IC à 95 % (~0.95) avec quantile de Student,
#      et sous-couverture si l'on utilise à tort le quantile normal (petit n) ;
#  (D) loi de la statistique t pivotale sous normalité = Student(n-p).
#
# DGP : y = X beta + eps, eps ~ N(0, sigma^2 I). X fixe (raisonnement
# conditionnel à X, cf. hypothèses H1-H5 de derivations/01_ols.qmd).
# =============================================================================

source("R/00_linalg.R")
source("R/01_ols.R")
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- DGP --------------------------------------------------------------------
n <- 15L                         # petit échantillon : t vs normal doit compter
p <- 4L
beta_true <- c(2, -1.5, 0, 0.8)  # 3e coef nul (test sous H0 vrai)
sigma <- 1.2
Xmat <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
colnames(Xmat) <- c("(Intercept)", "x1", "x2", "x3")
XtXinv <- chol2inv(chol(crossprod(Xmat)))
df_res <- n - p

# Estimateur alternatif linéaire sans biais : C = A + D avec D = D0 %*% M,
# DX = 0 (donc CX = I). Var(beta_tilde) = sigma^2 (X'X)^{-1} + sigma^2 DD'.
A  <- XtXinv %*% t(Xmat)                         # (X'X)^{-1} X'
M  <- diag(n) - Xmat %*% XtXinv %*% t(Xmat)      # I - H
set.seed(7)
D  <- matrix(rnorm(p * n), p, n) %*% M           # DX = 0
C_alt <- A + D

# ---- Boucle Monte Carlo -----------------------------------------------------
R <- 5000L
j_test <- 3L                     # coefficient suivi (x2, valeur vraie = 0)
z95 <- qnorm(0.975)
t95 <- qt(0.975, df_res)

beta_hat  <- matrix(NA_real_, R, p)
beta_alt  <- matrix(NA_real_, R, p)
cover_t   <- matrix(FALSE, R, p)   # IC de Student contient beta_vrai ?
cover_z   <- matrix(FALSE, R, p)   # IC "normal" (naïf) contient beta_vrai ?
t_pivot   <- numeric(R)

dat <- data.frame(Xmat[, -1]); names(dat) <- c("x1", "x2", "x3")

for (r in seq_len(R)) {
  y <- as.numeric(Xmat %*% beta_true + rnorm(n, sd = sigma))
  dat$y <- y
  fit <- ols_fit(y ~ x1 + x2 + x3, dat)
  b  <- fit$coefficients
  se <- sqrt(diag(fit$vcov))
  beta_hat[r, ] <- b
  beta_alt[r, ] <- as.numeric(C_alt %*% y)
  cover_t[r, ] <- (b - t95 * se <= beta_true) & (beta_true <= b + t95 * se)
  cover_z[r, ] <- (b - z95 * se <= beta_true) & (beta_true <= b + z95 * se)
  t_pivot[r] <- (b[j_test] - beta_true[j_test]) / se[j_test]   # ~ t_{n-p}
}

# ---- (A) Non-biais ----------------------------------------------------------
biais <- colMeans(beta_hat) - beta_true
cat("=== (A) Non-biais : moyenne(beta_hat) - beta_vrai ===\n")
print(round(setNames(biais, colnames(Xmat)), 4))

# ---- (B) BLUE : variances comparées -----------------------------------------
var_ols <- apply(beta_hat, 2, var)
var_alt <- apply(beta_alt, 2, var)
tabB <- data.frame(coef = colnames(Xmat),
                   var_OLS = round(var_ols, 4),
                   var_alt = round(var_alt, 4),
                   ratio = round(var_alt / var_ols, 3))
cat("\n=== (B) Gauss-Markov : Var(alternatif) >= Var(OLS) ===\n")
print(tabB, row.names = FALSE)

# ---- (C) Couverture ---------------------------------------------------------
cov_t <- colMeans(cover_t); cov_z <- colMeans(cover_z)
tabC <- data.frame(coef = colnames(Xmat),
                   couv_Student = round(cov_t, 3),
                   couv_normale = round(cov_z, 3))
cat("\n=== (C) Couverture empirique des IC à 95 % (R =", R, ") ===\n")
print(tabC, row.names = FALSE)
cat("Student -> ~0.95 (correct) ; normale -> sous-couvre (petit n =", n, ").\n")

# ---- (D) Loi de t : KS contre Student ---------------------------------------
ks <- suppressWarnings(ks.test(t_pivot, "pt", df = df_res))
cat("\n=== (D) t pivotal vs Student(", df_res, ") — test KS ===\n", sep = "")
cat("D =", round(ks$statistic, 4), " p-value =", round(ks$p.value, 3), "\n")

# ---- Graphiques -------------------------------------------------------------
gg1 <- ggplot(data.frame(t = t_pivot), aes(t)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 fill = "grey80", colour = "white") +
  stat_function(fun = dt, args = list(df = df_res), aes(colour = "Student(n-p)"),
                linewidth = 1) +
  stat_function(fun = dnorm, aes(colour = "Normale(0,1)"),
                linewidth = 1, linetype = "dashed") +
  scale_colour_manual(values = c("Student(n-p)" = "#1b7837", "Normale(0,1)" = "#d73027")) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(title = "Loi empirique de la statistique t sous normalité",
       subtitle = paste0("n = ", n, ", p = ", p,
                         " : queues plus épaisses que la normale (d'où le quantile de Student)"),
       x = expression((hat(beta)[j] - beta[j]) / hat(se)), y = "densité", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_01_t_student.png"), gg1, width = 8, height = 5, dpi = 120)

covdf <- rbind(data.frame(coef = colnames(Xmat), couv = cov_t, methode = "Student (correct)"),
               data.frame(coef = colnames(Xmat), couv = cov_z, methode = "Normale (naïf)"))
gg2 <- ggplot(covdf, aes(coef, couv, fill = methode)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.85, 1)) +
  labs(title = "Couverture empirique des IC à 95 %",
       subtitle = "Quantile de Student : ~0.95 ; quantile normal : sous-couvre en petit échantillon",
       x = NULL, y = "couverture", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_01_couverture.png"), gg2, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_01_t_student.png, mc_01_couverture.png\n", sep = "")

# ---- Interprétation ---------------------------------------------------------
# (A) biais ~ 0 : conforme au caractère sans biais (H1-H3).
# (B) ratio var_alt/var_OLS >= 1 pour chaque coef : illustration de Gauss-Markov
#     (OLS est le plus efficace parmi les linéaires sans biais).
# (C) la couverture Student ~0.95 ; la couverture normale < 0.95 car en petit
#     échantillon la vraie loi est Student, à queues plus épaisses.
# (D) KS ne rejette pas l'adéquation à Student(n-p) : la loi exacte est atteinte.
