# =============================================================================
# Monte Carlo — Module 41 : sur des formes NON CONVEXES, le k-means echoue ; le
# spectral et DBSCAN reussissent. On compare l'ARI (accord avec la verite) sur
# deux cercles concentriques et deux "lunes".
# =============================================================================

for (f in c("41_clustering", "11_pca_kmeans_em", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ari <- function(a, b) { tab <- table(a, b); s <- sum(choose(tab, 2))
  a1 <- sum(choose(rowSums(tab), 2)); b1 <- sum(choose(colSums(tab), 2))
  n <- length(a); e <- a1 * b1 / choose(n, 2); (s - e) / ((a1 + b1) / 2 - e) }

gen_circles <- function(n) { th <- runif(n, 0, 2*pi); r <- c(rep(1, n/2), rep(4, n/2)) + rnorm(n, 0, 0.12)
  list(X = cbind(r*cos(th), r*sin(th)), y = rep(1:2, each = n/2)) }
gen_moons <- function(n) { t1 <- runif(n/2, 0, pi); t2 <- runif(n/2, 0, pi)
  X <- rbind(cbind(cos(t1), sin(t1)), cbind(1 - cos(t2), 0.5 - sin(t2))) + matrix(rnorm(n*2, 0, 0.1), n, 2)
  list(X = X, y = rep(1:2, each = n/2)) }

R <- 40; A <- matrix(0, 2, 3, dimnames = list(c("cercles","lunes"), c("k-means","spectral","DBSCAN")))
for (r in seq_len(R)) {
  for (g in 1:2) {
    d <- if (g == 1) gen_circles(200) else gen_moons(200)
    A[g, 1] <- A[g, 1] + ari(kmeans(d$X, 2, nstart = 5)$cluster, d$y)
    A[g, 2] <- A[g, 2] + ari(spectral_clustering(d$X, 2, gamma = if (g==1) 2 else 8), d$y)
    db <- dbscan_fit(d$X, eps = if (g==1) 0.45 else 0.28, minPts = 4)
    A[g, 3] <- A[g, 3] + ari(db, d$y)
  }
}
A <- A / R
cat("=== Accord avec la verite (ARI) sur formes non convexes ===\n\n")
cat(sprintf("%-10s %10s %10s %10s\n", "forme", "k-means", "spectral", "DBSCAN"))
for (g in 1:2) cat(sprintf("%-10s %10.3f %10.3f %10.3f\n", rownames(A)[g], A[g,1], A[g,2], A[g,3]))
cat("\n=> Le k-means (groupes convexes) s'effondre sur cercles et lunes ; le\n")
cat("   spectral et DBSCAN, fondes sur la connectivite/densite, les recuperent.\n")

df <- data.frame(forme = rep(rownames(A), 3), methode = rep(colnames(A), each = 2), ARI = as.numeric(A))
gg <- ggplot(df, aes(methode, ARI, fill = methode)) + geom_col(show.legend = FALSE) +
  facet_wrap(~forme) + coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Clustering de formes non convexes : k-means vs spectral vs DBSCAN",
       subtitle = "ARI (1 = partition parfaite) ; le k-means echoue, les autres reussissent",
       x = NULL, y = "ARI") + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out_dir, "mc_41_clustering.png"), gg, width = 8.5, height = 5, dpi = 120)
cat("\nGraphique -> mc_41_clustering.png\n")
