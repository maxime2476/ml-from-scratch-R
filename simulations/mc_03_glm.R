# =============================================================================
# Monte Carlo — Module 3 : GLM logistique
#  (A) Normalité asymptotique du MLE : (beta_hat - beta)/se -> N(0,1) quand n
#      croît (la loi exacte n'est pas normale en petit n).
#  (B) Comparaison Wald / LR / score en petit échantillon : taille empirique du
#      test de H0: beta_j = 0 au niveau nominal 5 %, et puissance.
#
# DGP : y_i ~ Bernoulli(sigma(eta_i)), eta = b0 + b1 x1 + b2 x2.
# =============================================================================

for (f in c("00_linalg", "01_ols", "03_glm_irls"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

beta_true <- c(b0 = -0.4, b1 = 1.0, b2 = 0)   # b2 = 0 : coef nul (taille du test)

sim_one <- function(n) {
  x1 <- rnorm(n); x2 <- rnorm(n)
  eta <- beta_true[1] + beta_true[2] * x1 + beta_true[3] * x2
  y <- rbinom(n, 1, plogis(eta))
  data.frame(y = y, x1 = x1, x2 = x2)
}

# =============================================================================
# (A) Normalité asymptotique de beta_hat (coef b1)
# =============================================================================
ns <- c(25L, 50L, 200L)
R_A <- 3000L
tstats <- list()
for (n in ns) {
  tt <- numeric(0)
  reps <- 0L
  while (reps < R_A) {
    d <- sim_one(n)
    fit <- tryCatch(glm_irls(y ~ x1 + x2, d, "binomial"), warning = function(w) NULL,
                    error = function(e) NULL)
    if (is.null(fit) || any(!is.finite(fit$se))) next    # rejet des cas séparés
    reps <- reps + 1L
    tt <- c(tt, (fit$coefficients["x1"] - beta_true[2]) / fit$se["x1"])
  }
  tstats[[as.character(n)]] <- tt
}

skew   <- function(t) mean((t - mean(t))^3) / sd(t)^3
exkurt <- function(t) mean((t - mean(t))^4) / sd(t)^4 - 3

cat("=== (A) Normalité asymptotique : statistique (b1_hat - b1)/se ===\n")
tabA <- data.frame(
  n = ns,
  moyenne = sapply(tstats, mean),        # biais O(1/n) du MLE : -> 0 lentement
  ecart_type = sapply(tstats, sd),        # -> 1
  skewness = sapply(tstats, skew),        # -> 0 (forme symétrique)
  ex_kurtosis = sapply(tstats, exkurt))   # -> 0 (queues normales)
print(round(tabA, 3), row.names = FALSE)
cat("Forme -> N(0,1) : sd, skewness et kurtosis convergent. La moyenne tend vers\n",
    "0 plus lentement : c'est le biais O(1/n) du MLE logistique (petit échantillon).\n")

dfA <- do.call(rbind, lapply(names(tstats), function(nm)
  data.frame(t = tstats[[nm]], n = factor(paste0("n = ", nm), levels = paste0("n = ", ns)))))
ggA <- ggplot(dfA, aes(t)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "grey80", colour = "white") +
  stat_function(fun = dnorm, colour = "#d73027", linewidth = 1) +
  facet_wrap(~ n) + coord_cartesian(xlim = c(-4, 4)) +
  labs(title = "Normalité asymptotique du MLE logistique",
       subtitle = "Loi empirique de (b1_hat - b1)/se vs N(0,1) (rouge)",
       x = "statistique standardisée", y = "densité") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_03_normalite.png"), ggA, width = 9, height = 4, dpi = 120)

# =============================================================================
# (B) Wald / LR / score : taille (H0: b2 = 0) et puissance (b2 != 0)
# =============================================================================
R_B <- 3000L
Rmat <- matrix(c(0, 0, 1), 1, 3)   # b2 = 0
alpha <- 0.05; crit <- qchisq(1 - alpha, df = 1)

run_size_power <- function(n, b2) {
  bt <- c(beta_true[1:2], b2)
  rej <- c(Wald = 0, LR = 0, Score = 0); reps <- 0L
  while (reps < R_B) {
    x1 <- rnorm(n); x2 <- rnorm(n)
    y <- rbinom(n, 1, plogis(bt[1] + bt[2] * x1 + bt[3] * x2))
    d <- data.frame(y = y, x1 = x1, x2 = x2)
    full <- tryCatch(glm_irls(y ~ x1 + x2, d, "binomial"), warning = function(w) NULL, error = function(e) NULL)
    red  <- tryCatch(glm_irls(y ~ x1, d, "binomial"), warning = function(w) NULL, error = function(e) NULL)
    if (is.null(full) || is.null(red) || any(!is.finite(full$se))) next
    reps <- reps + 1L
    rej["Wald"]  <- rej["Wald"]  + (wald_test(full, Rmat)$statistic  > crit)
    rej["LR"]    <- rej["LR"]    + (lr_test(full, red)$statistic     > crit)
    rej["Score"] <- rej["Score"] + (score_test(full, red)$statistic  > crit)
  }
  rej / R_B
}

cat("\n=== (B) Taille empirique du test de H0: b2 = 0 (nominal 5 %) ===\n")
sizeTab <- data.frame(n = integer(0), Wald = numeric(0), LR = numeric(0), Score = numeric(0))
for (n in c(20L, 40L, 100L)) {
  s <- run_size_power(n, b2 = 0)
  sizeTab <- rbind(sizeTab, data.frame(n = n, Wald = s["Wald"], LR = s["LR"], Score = s["Score"]))
}
print(round(sizeTab, 3), row.names = FALSE)
cat("En petit n, les trois tailles diffèrent (Wald souvent distordu) ; elles\n",
    "convergent vers 5 % quand n croît (équivalence asymptotique, Th. 3.3).\n")

sizeLong <- reshape(sizeTab, varying = c("Wald", "LR", "Score"), v.names = "taille",
                    timevar = "test", times = c("Wald", "LR", "Score"), direction = "long")
ggB <- ggplot(sizeLong, aes(factor(n), taille, fill = test)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(title = "Taille empirique des tests Wald / LR / score (H0: b2 = 0)",
       subtitle = "Nominal 5 % (pointillé) ; distorsions en petit n, convergence quand n croît",
       x = "n", y = "taux de rejet sous H0", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_03_wald_lr_score.png"), ggB, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_03_normalite.png, mc_03_wald_lr_score.png\n", sep = "")
