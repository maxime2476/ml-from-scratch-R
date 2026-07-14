# =============================================================================
# Monte Carlo — Module 45 : le GARCH capture le REGROUPEMENT de volatilite.
# (1) La volatilite conditionnelle estimee suit la volatilite realisee.
# (2) Ignorer l'ARCH sous-estime le risque : la VaR d'un modele homoscedastique
#     est trop optimiste en periode agitee.
# =============================================================================

for (f in c("45_garch", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sim_garch <- function(n, w = 0.05, al = 0.12, be = 0.85) {
  x <- numeric(n); s2 <- numeric(n); s2[1] <- w / (1 - al - be)
  for (t in 2:n) { s2[t] <- w + al * x[t - 1]^2 + be * s2[t - 1]; x[t] <- rnorm(1) * sqrt(s2[t]) }
  list(x = x, sigma = sqrt(s2))
}

## (1) suivi de la volatilite ------------------------------------------------
set.seed(1); d <- sim_garch(1500); m <- garch_fit(d$x)
cat(sprintf("=== (1) GARCH(1,1) estime (vrais : omega=0.05, alpha=0.12, beta=0.85) ===\n"))
cat(sprintf("  omega %.3f | alpha %.3f | beta %.3f | persistance %.3f\n", m$omega, m$alpha, m$beta, m$persistence))
cat(sprintf("  correlation sigma estime vs sigma VRAI : %.3f\n\n", cor(m$sigma, d$sigma)))

## (2) couverture d'une VaR a 1 % : GARCH vs modele homoscedastique ----------
R <- 300; hit_h <- hit_g <- 0; ntot <- 0; z <- qnorm(0.01)
for (r in seq_len(R)) {
  s <- sim_garch(600); x <- s$x; te <- 301:600
  # VaR 1% (borne inferieure) : homoscedastique (sd constant) vs GARCH (sigma_t)
  sd_h <- sd(x[1:300])
  mg <- garch_fit(x[1:300])                                # refit sur la 1re moitie -> sigma persistant
  sig_g <- mg$sigma[300]                                   # derniere volatilite (approx pour la suite)
  for (t in te) {
    hit_h <- hit_h + (x[t] < z * sd_h)                     # depassement de la VaR homoscedastique
    hit_g <- hit_g + (x[t] < z * sqrt(mg$omega + mg$alpha * x[t-1]^2 + mg$beta * sig_g^2))
    sig_g <- sqrt(mg$omega + mg$alpha * x[t-1]^2 + mg$beta * sig_g^2)
  }
  ntot <- ntot + length(te)
}
cat("=== (2) Depassements d'une VaR a 1 % (nominal 0.01) ===\n")
cat(sprintf("  VaR homoscedastique (sd constant) : %.3f  (mal calibree)\n", hit_h / ntot))
cat(sprintf("  VaR GARCH (volatilite dynamique)  : %.3f  (~ nominal)\n", hit_g / ntot))
cat("\n=> Le GARCH suit la volatilite (regroupement) ; ignorer l'ARCH -> une VaR\n")
cat("   mal calibree (trop de depassements quand la volatilite monte).\n")

df <- data.frame(t = seq_along(d$x), sig_vrai = d$sigma, sig_est = m$sigma)
gg <- ggplot(df, aes(t)) +
  geom_line(aes(y = sig_vrai, colour = "vraie"), linewidth = 0.5) +
  geom_line(aes(y = sig_est, colour = "GARCH estimee"), linewidth = 0.5) +
  labs(title = "Volatilite conditionnelle : le GARCH suit le regroupement",
       subtitle = "volatilite vraie vs estimee (sigma_t) le long de la serie",
       x = "temps", y = "volatilite conditionnelle", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_45_garch.png"), gg, width = 9, height = 5, dpi = 120)
cat("\nGraphique -> mc_45_garch.png\n")
