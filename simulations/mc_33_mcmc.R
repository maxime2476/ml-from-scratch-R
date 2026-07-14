# =============================================================================
# Monte Carlo — Module 33 : le pas de proposition gouverne le melange
# Un pas trop PETIT : acceptation ~1 mais chaine lente (autocorrelation elevee,
# ESS faible). Un pas trop GRAND : rejets frequents, ESS faible aussi. Il existe
# un pas OPTIMAL (taux d'acceptation ~ 0.2-0.4) maximisant l'ESS.
# =============================================================================

for (f in c("33_mcmc", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# cible : normale standard multivariee (dim 1 ici)
lt <- function(x) dnorm(x, 0, 1, log = TRUE)
psds <- c(0.1, 0.3, 0.7, 1.5, 3, 6, 12); R <- 40; N <- 8000
res <- matrix(0, length(psds), 2)                          # acceptation, ESS
for (r in seq_len(R)) for (j in seq_along(psds)) {
  m <- metropolis_hastings(lt, 0, psds[j], N)
  ch <- m$chain[2001:N, 1]
  res[j, 1] <- res[j, 1] + m$accept_rate
  res[j, 2] <- res[j, 2] + ess(ch)
}
res <- res / R
cat("=== Metropolis-Hastings : pas de proposition vs melange (cible N(0,1)) ===\n\n")
cat(sprintf("%8s %12s %10s\n", "pas", "acceptation", "ESS"))
for (j in seq_along(psds)) cat(sprintf("%8.1f %12.2f %10.0f\n", psds[j], res[j, 1], res[j, 2]))
jbest <- which.max(res[, 2])
cat(sprintf("\n=> ESS maximale au pas %.1f (acceptation %.2f) : ni trop petit (chaine\n",
            psds[jbest], res[jbest, 1]))
cat("   lente) ni trop grand (rejets). Le taux d'acceptation optimal est modere.\n")

df <- data.frame(pas = psds, acceptation = res[, 1], ESS = res[, 2])
gg <- ggplot(df, aes(acceptation, ESS)) +
  geom_line(colour = "#00798c") + geom_point(aes(size = pas), colour = "#d1495b") +
  labs(title = "MCMC : la taille effective culmine a un taux d'acceptation modere",
       subtitle = "trop de rejets (pas grand) ou trop d'autocorrelation (pas petit) tuent l'ESS",
       x = "taux d'acceptation", y = "taille d'echantillon effective (ESS)", size = "pas") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_33_mcmc.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_33_mcmc.png\n")
