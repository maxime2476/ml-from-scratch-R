# =============================================================================
# Monte Carlo — Module 21 : panel / effets fixes
#  (A) Sous effet fixe alpha_i CORRÉLÉ à x, l'OLS groupé est biaisé et sous-couvre ;
#      l'estimateur within (effets fixes) est sans biais et couvre.
#  (B) Sous corrélation SÉRIELLE (AR1), les SE within classiques sous-couvrent ;
#      les SE groupées (clustered) restaurent la couverture.
# =============================================================================

for (f in c("00_linalg", "01_ols", "21_panel", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

b_true <- 2

# =============================================================================
# (A) Biais : OLS groupé vs effets fixes
# =============================================================================
N <- 100; T <- 5; R <- 2000
b_pool <- b_fe <- numeric(R); cov_pool <- cov_fe <- logical(R)
for (r in seq_len(R)) {
  alpha <- rnorm(N)
  id <- rep(seq_len(N), each = T)
  x <- 0.8 * alpha[id] + rnorm(N * T)            # x corrélé à l'effet fixe
  y <- 1 + b_true * x + alpha[id] + rnorm(N * T)
  d <- data.frame(id = id, x = x, y = y)
  ml <- lm(y ~ x, d); b_pool[r] <- coef(ml)["x"]
  se_p <- summary(ml)$coefficients["x", "Std. Error"]
  cov_pool[r] <- abs(b_pool[r] - b_true) <= 1.96 * se_p
  fe <- fe_within(y ~ x, d, "id"); b_fe[r] <- fe$coefficients["x"]
  cov_fe[r] <- abs(b_fe[r] - b_true) <= 1.96 * fe$se["x"]
}
sp <- mc_summary(b_pool, b_true); sf <- mc_summary(b_fe, b_true)
cat("=== (A) OLS groupé vs effets fixes (vrai =", b_true, ") ===\n")
cat(sprintf("OLS groupé   : biais %+.3f (%.3f)  couverture %.3f\n", sp$bias, sp$bias_se, mean(cov_pool)))
cat(sprintf("Effets fixes : biais %+.3f (%.3f)  couverture %.3f\n", sf$bias, sf$bias_se, mean(cov_fe)))

dfA <- rbind(data.frame(b = b_pool, est = "OLS groupé"), data.frame(b = b_fe, est = "effets fixes"))
ggA <- ggplot(dfA, aes(b, fill = est)) + geom_density(alpha = 0.5) +
  geom_vline(xintercept = b_true, linetype = "dashed") +
  labs(title = "Panel : OLS groupé (biaisé) vs effets fixes (sans biais)",
       subtitle = "L'effet fixe alpha_i corrélé à x biaise l'OLS groupé ; le within l'élimine",
       x = expression(hat(beta)), y = "densité", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_21_biais.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) SE naïves vs groupées sous corrélation sérielle (AR1)
# =============================================================================
N2 <- 80; T2 <- 10; rho <- 0.8; R2 <- 2000
ar1 <- function(T) { e <- numeric(T); e[1] <- rnorm(1); for (t in 2:T) e[t] <- rho * e[t-1] + rnorm(1); e }
cov_naive <- cov_clust <- logical(R2)
for (r in seq_len(R2)) {
  id <- rep(seq_len(N2), each = T2)
  x <- as.numeric(sapply(seq_len(N2), function(i) ar1(T2)))
  eps <- as.numeric(sapply(seq_len(N2), function(i) ar1(T2)))
  y <- 1 + b_true * x + rnorm(N2)[id] + eps
  d <- data.frame(id = id, x = x, y = y)
  fe <- fe_within(y ~ x, d, "id")
  cov_naive[r] <- abs(fe$coefficients["x"] - b_true) <= 1.96 * fe$se["x"]
  cov_clust[r] <- abs(fe$coefficients["x"] - b_true) <= 1.96 * fe$se_cluster["x"]
}
cn <- coverage_mc(cov_naive); cc <- coverage_mc(cov_clust)
cat("\n=== (B) Couverture des IC à 95 % sous corrélation sérielle (AR1, rho=0.8) ===\n")
cat(sprintf("SE within classiques : %.3f ± %.3f  (SOUS-COUVRE)\n", cn$coverage, cn$se))
cat(sprintf("SE groupées (cluster): %.3f ± %.3f  (restaure ~0.95)\n", cc$coverage, cc$se))
cat("\nGraphique -> ", file.path(out_dir, "mc_21_biais.png"), "\n")
