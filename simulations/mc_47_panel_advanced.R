# =============================================================================
# Monte Carlo — Module 47 : (1) le biais de Nickell des effets fixes dynamiques
# et sa correction par IV ; (2) le controle synthetique reproduit le contrefactuel.
# =============================================================================

for (f in c("47_panel_advanced", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) Biais de Nickell vs IV, selon la longueur du panel T -------------------
rho <- 0.5; N <- 150; Ts <- c(4, 6, 10, 20); R <- 60
tab <- matrix(0, length(Ts), 2, dimnames = list(paste0("T=", Ts), c("within (Nickell)", "IV")))
for (k in seq_along(Ts)) {
  T <- Ts[k]
  for (r in seq_len(R)) {
    set.seed(1000 * k + r); a <- rnorm(N); Y <- matrix(0, T, N)
    for (i in 1:N) { yi <- rnorm(1); for (t in 1:T) { yi <- rho * yi + a[i] + rnorm(1); Y[t, i] <- yi } }
    d <- data.frame(id = rep(1:N, each = T), time = rep(1:T, N), y = as.numeric(Y))
    tab[k, 1] <- tab[k, 1] + dynamic_panel_fe(d)
    tab[k, 2] <- tab[k, 2] + dynamic_panel_iv(d)$rho
  }
}
tab <- tab / R
cat(sprintf("=== (1) Estimation de rho = %.1f (panel dynamique) ===\n\n", rho))
cat(sprintf("%8s %18s %10s\n", "T", "within (Nickell)", "IV"))
for (k in seq_along(Ts)) cat(sprintf("%8d %18.3f %10.3f\n", Ts[k], tab[k, 1], tab[k, 2]))
cat("\n=> Le within est BIAISE vers le bas (biais de Nickell ~ -1/T, severe en panel\n")
cat("   court) ; l'IV (Anderson-Hsiao) le corrige a toute longueur.\n\n")

## (2) Controle synthetique : contrefactuel + effet ---------------------------
set.seed(1); Tt <- 40; J <- 12; pre <- 1:30; post <- 31:40
Y0 <- matrix(0, Tt, J); for (j in 1:J) Y0[, j] <- cumsum(rnorm(Tt, 0.05)) + rnorm(Tt, 0, 0.3)
w <- c(0.4, 0.3, 0.2, 0.1, rep(0, J - 4)); Y1 <- as.numeric(Y0 %*% w) + rnorm(Tt, 0, 0.1)
effet <- 3; Y1[post] <- Y1[post] + effet                     # intervention en t=31
sc <- synthetic_control(Y1, Y0, pre)
cat("=== (2) Controle synthetique (effet injecte = +3 en post) ===\n")
cat(sprintf("  ajustement pre-traitement (RMSE) : %.3f\n", sqrt(mean((Y1[pre] - sc$synthetic[pre])^2))))
cat(sprintf("  effet moyen post-traitement estime : %.3f  (vrai +%d)\n", mean(sc$effect), effet))
cat("=> Le synthetique colle avant l'intervention puis diverge de l'effet injecte.\n")

# figures
df1 <- data.frame(T = rep(Ts, 2), rho = c(tab[,1], tab[,2]),
                  methode = rep(c("within (Nickell)", "IV (Anderson-Hsiao)"), each = length(Ts)))
g1 <- ggplot(df1, aes(T, rho, colour = methode)) + geom_line(linewidth = 1) + geom_point() +
  geom_hline(yintercept = rho, linetype = "dashed") +
  labs(title = "Panel dynamique : biais de Nickell (within) vs correction IV",
       subtitle = paste0("vrai rho = ", rho, " (tirete) ; le within converge vers un mauvais point"),
       x = "longueur du panel T", y = "rho estime", colour = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_47_nickell.png"), g1, width = 8, height = 5, dpi = 120)

df2 <- data.frame(t = 1:Tt, traite = Y1, synthetique = sc$synthetic)
g2 <- ggplot(df2, aes(t)) + geom_line(aes(y = traite, colour = "unite traitee"), linewidth = 1) +
  geom_line(aes(y = synthetique, colour = "controle synthetique"), linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 30.5, linetype = "dotted") +
  labs(title = "Controle synthetique : le contrefactuel colle avant, diverge apres",
       subtitle = "l'ecart post-intervention (t > 30) estime l'effet causal",
       x = "temps", y = "resultat", colour = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_47_synth.png"), g2, width = 8, height = 5, dpi = 120)
cat("\nGraphiques -> mc_47_nickell.png, mc_47_synth.png\n")
