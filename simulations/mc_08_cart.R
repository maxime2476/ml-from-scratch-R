# =============================================================================
# Monte Carlo — Module 8 : CART
#  (A) Gini vs entropie : les premiers splits sont quasi identiques (Prop. 8.1).
#  (B) Contre-exemple XOR : le greedy est myope (profondeur 1 inutile), résolu
#      en profondeur 2.
#  (C) Erreur vs complexité : sur/sous-ajustement selon la profondeur et
#      l'élagage coût-complexité (biais-variance de l'arbre).
# =============================================================================

for (f in c("00_linalg", "06_validation", "08_cart"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# (A) Gini vs entropie : premiers splits
# =============================================================================
R <- 200; agree_var <- 0; thr_diff <- numeric(R)
for (r in seq_len(R)) {
  n <- 300
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- factor(ifelse(0.9 * d$x1 - 0.6 * d$x2 + rnorm(n) > 0, "A", "B"))
  fg <- cart_fit(y ~ x1 + x2 + x3, d, "class", kind = "gini", min_split = 20)
  fe <- cart_fit(y ~ x1 + x2 + x3, d, "class", kind = "entropy", min_split = 20)
  agree_var <- agree_var + (fg$tree$var == fe$tree$var)
  thr_diff[r] <- if (fg$tree$var == fe$tree$var) abs(fg$tree$val - fe$tree$val) else NA
}
cat("=== (A) Gini vs entropie (premiers splits, R =", R, ") ===\n")
cat(sprintf("Même variable de split : %.1f %% des cas\n", 100 * agree_var / R))
cat(sprintf("Écart de seuil quand même variable : médiane %.4f\n", median(thr_diff, na.rm = TRUE)))
cat("=> Gini et entropie produisent des splits quasi identiques (Prop. 8.1).\n")

# =============================================================================
# (B) Contre-exemple XOR
# =============================================================================
m <- 150
g <- expand.grid(x1 = c(0, 1), x2 = c(0, 1))
dxor <- g[rep(seq_len(4), each = m), ]
dxor$y <- factor(as.integer(xor(dxor$x1 == 1, dxor$x2 == 1)))
accXOR <- sapply(0:3, function(dep) {
  fit <- cart_fit(y ~ x1 + x2, dxor, "class", max_depth = dep,
                  min_split = 2, min_leaf = 1, min_gain = -1)  # accepte splits à gain nul
  mean(predict_cart(fit, dxor) == dxor$y)
})
cat("\n=== (B) XOR : précision selon la profondeur maximale ===\n")
print(data.frame(profondeur = 0:3, precision = round(accXOR, 3)), row.names = FALSE)
cat("Profondeur 0-1 : 50 % (greedy myope) ; profondeur >= 2 : 100 % (structure révélée).\n")

# =============================================================================
# (C) Erreur vs complexité (régression) : profondeur et élagage
# =============================================================================
f_true <- function(X) sin(2 * X[, 1]) + 0.5 * X[, 2]^2
gen <- function(n) { X <- matrix(runif(n * 2, -2, 2), n, 2)
                     list(X = X, y = f_true(X) + rnorm(n, sd = 0.4)) }
tr <- gen(300); te <- gen(3000)
dtr <- data.frame(x1 = tr$X[, 1], x2 = tr$X[, 2], y = tr$y)
dte <- data.frame(x1 = te$X[, 1], x2 = te$X[, 2], y = te$y)

depths <- 1:12
errC <- data.frame(depth = depths, train = NA, test = NA)
for (i in seq_along(depths)) {
  fit <- cart_fit(y ~ x1 + x2, dtr, "anova", max_depth = depths[i], min_split = 5, min_leaf = 2)
  errC$train[i] <- mean((predict_cart(fit, dtr) - dtr$y)^2)
  errC$test[i]  <- mean((predict_cart(fit, dte) - dte$y)^2)
}
cat("\n=== (C) Erreur train/test selon la profondeur (régression) ===\n")
print(round(errC, 4), row.names = FALSE)
d_opt <- depths[which.min(errC$test)]
cat(sprintf("Profondeur* (test min) = %d ; au-delà, sur-ajustement (train baisse, test remonte).\n", d_opt))

# Élagage coût-complexité : erreur test selon alpha (arbre profond élagué)
big <- cart_fit(y ~ x1 + x2, dtr, "anova", max_depth = 20, min_split = 5, min_leaf = 2)
alphas <- 10^seq(-4, -0.5, length.out = 15)
pruneC <- data.frame(alpha = alphas, leaves = NA, test = NA)
for (i in seq_along(alphas)) {
  pr <- cost_complexity_prune(big, alphas[i])
  pruneC$leaves[i] <- n_leaves(pr)
  pruneC$test[i]   <- mean((predict_cart(pr, dte) - dte$y)^2)
}
cat("\n=== Élagage coût-complexité : erreur test selon alpha ===\n")
print(round(pruneC, 4), row.names = FALSE)

gg1 <- ggplot(reshape(errC, varying = c("train", "test"), v.names = "err",
                      timevar = "ensemble", times = c("train", "test"), direction = "long"),
              aes(depth, err, colour = ensemble)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  geom_vline(xintercept = d_opt, linetype = "dashed", colour = "grey40") +
  labs(title = "CART : erreur train/test selon la profondeur",
       subtitle = "Le train baisse toujours ; le test est en U (sur-ajustement en profondeur)",
       x = "profondeur maximale", y = "EQM", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_08_profondeur.png"), gg1, width = 8, height = 5, dpi = 120)

gg2 <- ggplot(pruneC, aes(leaves, test)) +
  geom_line(linewidth = 1, colour = "#2166ac") + geom_point(size = 1.8, colour = "#2166ac") +
  labs(title = "CART : erreur test selon la taille de l'arbre élagué",
       subtitle = "Cost-complexity pruning : compromis taille (biais) / erreur (variance)",
       x = "nombre de feuilles", y = "EQM de test") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_08_elagage.png"), gg2, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir, "/mc_08_profondeur.png, mc_08_elagage.png\n", sep = "")
