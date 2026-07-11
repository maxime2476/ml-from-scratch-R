# =============================================================================
# Monte Carlo / illustrations — Module 11 : non supervisé
#  (A) ACP : les deux voies (variance/Lagrange et SVD) coïncident exactement.
#  (B) k-means = limite de l'EM : quand sigma² -> 0, les responsabilités
#      deviennent dures et la partition de l'EM sphérique = celle du k-means.
#  (C) Sensibilité à l'initialisation : k-means et EM ont des optima locaux.
# =============================================================================

for (f in c("00_linalg", "11_pca_kmeans_em")) source(file.path("R", paste0(f, ".R")))
suppressMessages({library(ggplot2); library(MASS)})

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) ACP : coïncidence des deux voies
# =============================================================================
Sigma_true <- matrix(c(4, 1.5, 1.5, 1), 2, 2)
X <- MASS::mvrnorm(500, c(0, 0), Sigma_true)
pc <- pca_fit(X)
eig <- eigen(cov(X))
cat("=== (A) ACP : voie variance (eigen de S) vs voie SVD ===\n")
cat("valeurs propres (voie 1) :", round(eig$values, 4), "\n")
cat("sdev^2       (voie 2 SVD) :", round(pc$sdev^2, 4), "\n")
cat("max|diff valeurs| :", max(abs(eig$values - pc$sdev^2)),
    " max|diff axes| :", max(abs(abs(eig$vectors) - abs(pc$rotation))), "\n")
cat("part de variance expliquée :", round(pc$var_explained, 3), "\n")

dfA <- as.data.frame(X); names(dfA) <- c("x1", "x2")
sc <- 2.5 * pc$sdev
axes <- data.frame(xend = sc * pc$rotation[1, ], yend = sc * pc$rotation[2, ])
ggA <- ggplot(dfA, aes(x1, x2)) + geom_point(alpha = 0.25, size = 1) +
  geom_segment(data = axes, aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#d73027", linewidth = 1) +
  coord_equal() +
  labs(title = "ACP : les deux voies donnent les mêmes axes principaux",
       subtitle = "Vecteurs propres de la covariance = vecteurs singuliers droits (Prop. 11.2)",
       x = expression(x[1]), y = expression(x[2])) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_11_pca.png"), ggA, width = 6, height = 5, dpi = 120)

# =============================================================================
# (B) k-means comme limite de l'EM sphérique (sigma^2 -> 0)
# =============================================================================
Xc <- rbind(MASS::mvrnorm(120, c(0, 0), diag(2)),
            MASS::mvrnorm(120, c(4, 4), diag(2)),
            MASS::mvrnorm(120, c(0, 5), diag(2)))
km <- kmeans_fit(Xc, 3, seed = 1)
centres <- km$centers

# Responsabilités sphériques (variance commune sigma^2 I) aux centres du k-means.
spherical_resp <- function(X, mu, sigma2) {
  logd <- sapply(seq_len(nrow(mu)), function(k)
    -colSums((t(X) - mu[k, ])^2) / (2 * sigma2))
  ls <- .logsumexp_rows(logd)
  exp(logd - ls)
}
sig2s <- 10^seq(1, -2, length.out = 12)
tabB <- data.frame(sigma2 = sig2s, resp_max = NA, accord_kmeans = NA)
for (i in seq_along(sig2s)) {
  g <- spherical_resp(Xc, centres, sig2s[i])
  tabB$resp_max[i] <- mean(apply(g, 1, max))
  tabB$accord_kmeans[i] <- mean(max.col(g, "first") == km$cluster)
}
cat("\n=== (B) k-means comme limite EM (sigma^2 -> 0) ===\n")
print(round(tabB, 4), row.names = FALSE)
cat("Quand sigma^2 -> 0 : responsabilité max -> 1 (affectation dure) et accord\n",
    "avec le k-means -> 100 %.\n", sep = "")

ggB <- ggplot(tabB, aes(sigma2)) +
  geom_line(aes(y = resp_max, colour = "responsabilité max moyenne"), linewidth = 1) +
  geom_point(aes(y = resp_max, colour = "responsabilité max moyenne")) +
  geom_line(aes(y = accord_kmeans, colour = "accord avec k-means"), linewidth = 1) +
  geom_point(aes(y = accord_kmeans, colour = "accord avec k-means")) +
  scale_x_log10() +
  labs(title = "k-means comme limite de l'EM gaussien sphérique",
       subtitle = "sigma^2 -> 0 : responsabilités dures et partition = k-means",
       x = expression(sigma^2~"(log)"), y = "proportion", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_11_kmeans_limite.png"), ggB, width = 8, height = 5, dpi = 120)

# =============================================================================
# (C) Sensibilité à l'initialisation
# =============================================================================
set.seed(9)
# Clusters proches et chevauchants + demande de K=6 (> vrai) : nombreux minima locaux.
Xs <- rbind(MASS::mvrnorm(80, c(0, 0), diag(2)),
            MASS::mvrnorm(80, c(2, 0), diag(2)),
            MASS::mvrnorm(80, c(0, 2), diag(2)),
            MASS::mvrnorm(80, c(2, 2), diag(2)),
            MASS::mvrnorm(80, c(1, 1), diag(2)))
Kc <- 6; nrun <- 200
inerties <- replicate(nrun, kmeans_fit(Xs, Kc, nstart = 1)$tot_withinss)
cat("\n=== (C) Sensibilité à l'initialisation (k-means, K=6, 200 départs) ===\n")
cat(sprintf("inertie : min = %.2f  médiane = %.2f  max = %.2f\n",
            min(inerties), median(inerties), max(inerties)))
cat(sprintf("proportion de départs atteignant le meilleur optimum (à 1%%) : %.1f %%\n",
            100 * mean(inerties <= 1.01 * min(inerties))))
cat("=> l'algorithme reste piégé dans des minima locaux : d'où les redémarrages.\n")

ggC <- ggplot(data.frame(inertie = inerties), aes(inertie)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "white") +
  geom_vline(xintercept = min(inerties), colour = "#1b7837", linetype = "dashed") +
  labs(title = "k-means : sensibilité à l'initialisation",
       subtitle = "200 départs aléatoires : plusieurs minima locaux (vert = meilleur trouvé)",
       x = "inertie intra-classe finale", y = "nombre de départs") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_11_init.png"), ggC, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_11_pca.png, mc_11_kmeans_limite.png, mc_11_init.png\n", sep = "")
