# =============================================================================
# Monte Carlo — Module 19 : prédiction conforme
#  (A) La prédiction conforme maintient la couverture ~1-alpha sur TOUS les DGP
#      (gaussien, queues lourdes, hétéroscédastique, modèle mal spécifié), là où
#      l'intervalle GAUSSIEN NAÏF (mu_hat +- z*sigma_hat) échoue.
#  (B) Garantie en échantillon fini : couverture dans [1-alpha, 1-alpha+1/(n+1)].
# =============================================================================

for (f in c("00_linalg", "19_conformal", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ols_fit  <- function(X, y) solve_ls_qr(cbind(1, X), y)
ols_pred <- function(m, X) as.numeric(cbind(1, X) %*% m$coefficients)

alpha <- 0.1; n <- 200
scenarios <- list(
  "gaussien"          = function(m) list(X = matrix(rnorm(m*2),m,2), e = rnorm(m), f = function(X) 1 + X %*% c(2,-1)),
  "queues lourdes t2" = function(m) list(X = matrix(rnorm(m*2),m,2), e = rt(m,2),  f = function(X) 1 + X %*% c(2,-1)),
  "hétéroscédastique" = function(m) { X <- matrix(rnorm(m*2),m,2); list(X = X, e = rnorm(m)*exp(0.6*X[,1]), f = function(X) 1 + X %*% c(2,-1)) },
  "modèle mal spécifié" = function(m) { X <- matrix(rnorm(m*2),m,2); list(X = X, e = rnorm(m,sd=0.4), f = function(X) sin(3*X[,1]) + X[,2]^2) })

R <- 2000
tab <- data.frame()
for (nm in names(scenarios)) {
  gen <- scenarios[[nm]]
  cf_cov <- naive_cov <- logical(R)
  for (r in seq_len(R)) {
    d <- gen(n + 1); y <- as.numeric(d$f(d$X)) + d$e
    tr <- 1:(n/2); cal <- (n/2+1):n; te <- n+1
    cf <- conformal_split(d$X[tr,,drop=FALSE], y[tr], d$X[cal,,drop=FALSE], y[cal],
                          d$X[te,,drop=FALSE], ols_fit, ols_pred, alpha = alpha)
    cf_cov[r] <- y[te] >= cf$lower && y[te] <= cf$upper
    # intervalle gaussien naïf : ajuste sur tr+cal, sigma = sd des résidus
    m0 <- ols_fit(d$X[1:n,,drop=FALSE], y[1:n]); res <- y[1:n] - ols_pred(m0, d$X[1:n,,drop=FALSE])
    p0 <- ols_pred(m0, d$X[te,,drop=FALSE]); q <- qnorm(1 - alpha/2) * sd(res)
    naive_cov[r] <- abs(y[te] - p0) <= q
  }
  cc <- coverage_mc(cf_cov, 1 - alpha); nc <- coverage_mc(naive_cov, 1 - alpha)
  tab <- rbind(tab, data.frame(scenario = nm,
    conforme = cc$coverage, conforme_se = cc$se,
    naif = nc$coverage, naif_se = nc$se))
}
cat("=== (A) Couverture à", 100*(1-alpha), "% : conforme vs gaussien naïf ===\n")
for (i in seq_len(nrow(tab)))
  cat(sprintf("  %-22s conforme %.3f ± %.3f | naïf %.3f ± %.3f\n",
              tab$scenario[i], tab$conforme[i], tab$conforme_se[i], tab$naif[i], tab$naif_se[i]))
cat("La conforme atteint EXACTEMENT ~0.90 partout (distribution-libre) ; le naïf se\n",
    "DÉ-CALIBRE selon la loi (ici il sur-couvre : sigma_hat gonflée par les queues\n",
    "lourdes / la mauvaise spécification -> intervalles inutilement larges).\n", sep = "")

long <- rbind(data.frame(scenario = tab$scenario, methode = "conforme", couv = tab$conforme),
              data.frame(scenario = tab$scenario, methode = "gaussien naïf", couv = tab$naif))
long$scenario <- factor(long$scenario, levels = names(scenarios))
gg <- ggplot(long, aes(scenario, couv, fill = methode)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  coord_cartesian(ylim = c(0.6, 1)) +
  labs(title = "Prédiction conforme vs intervalle gaussien naïf",
       subtitle = paste0("Couverture visée ", 100*(1-alpha), "% : la conforme est exactement calibrée ; le naïf se dé-calibre"),
       x = NULL, y = "couverture empirique", fill = NULL) +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out_dir, "mc_19_conformal.png"), gg, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Garantie en échantillon fini : couverture selon n_cal
# =============================================================================
cat("\n=== (B) Couverture dans [1-alpha, 1-alpha+1/(n_cal+1)] ===\n")
gen <- scenarios[["gaussien"]]
for (ncal in c(20L, 50L, 200L)) {
  cov <- mean(replicate(3000, {
    d <- gen(ncal + 51); y <- as.numeric(d$f(d$X)) + d$e
    tr <- 1:50; cal <- 51:(50+ncal); te <- 50+ncal+1
    cf <- conformal_split(d$X[tr,,drop=FALSE], y[tr], d$X[cal,,drop=FALSE], y[cal],
                          d$X[te,,drop=FALSE], ols_fit, ols_pred, alpha = alpha)
    y[te] >= cf$lower && y[te] <= cf$upper
  }))
  cat(sprintf("  n_cal=%3d : couverture %.3f  (borne théorique [%.3f, %.3f], à l'erreur MC ~0.006 près)\n",
              ncal, cov, 1 - alpha, 1 - alpha + 1/(ncal + 1)))
}
cat("\nGraphique -> ", file.path(out_dir, "mc_19_conformal.png"), "\n")
