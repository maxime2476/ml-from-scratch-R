# =============================================================================
# Monte Carlo / illustrations — Module 12 : MLP
#  (A) Courbe d'apprentissage et effet du taux d'apprentissage (perte vs époques).
#  (B) Effet de la largeur cachée d1 (sous-ajustement vs capacité).
#  (C) Précision du gradient numérique : erreur relative vs eps (O(eps²) puis
#      plancher d'arrondi) — illustration de la vérification du gradient.
# =============================================================================

source("R/12_mlp.R")
suppressMessages(library(ggplot2))

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# DGP de régression non linéaire
f_true <- function(x) sin(2 * x) + 0.3 * x
gen <- function(n, seed) { set.seed(seed)
  x <- matrix(runif(n, -3, 3), n, 1); list(x = x, y = f_true(x) + rnorm(n, sd = 0.15)) }
tr <- gen(300, 1); te <- gen(1000, 2)

# =============================================================================
# (A) Courbe d'apprentissage selon le taux d'apprentissage
# =============================================================================
lrs <- c(0.005, 0.02, 0.1, 0.4)
E <- 400
dfA <- do.call(rbind, lapply(lrs, function(lr) {
  m <- mlp_fit(tr$x, tr$y, hidden = 20, activation = "tanh", loss = "mse",
               epochs = E, lr = lr, batch = 32, seed = 1)
  data.frame(epoch = 1:E, loss = m$loss_hist, lr = factor(lr))
}))
cat("=== (A) Effet du taux d'apprentissage (perte train finale) ===\n")
for (lr in lrs) {
  h <- dfA$loss[dfA$lr == as.character(lr)]
  cat(sprintf("lr=%-5s : perte finale = %.4f %s\n", as.character(lr), tail(h, 1),
              if (any(!is.finite(h)) || tail(h, 1) > h[1]) "(instable)" else ""))
}
cat("Petit lr : descente lente ; grand lr : rapide mais peut osciller/diverger.\n")

ggA <- ggplot(dfA, aes(epoch, loss, colour = lr)) + geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(title = "MLP : courbe d'apprentissage selon le taux d'apprentissage",
       subtitle = "Petit lr converge lentement ; grand lr converge vite mais peut osciller",
       x = "époque", y = "perte d'entraînement (log)", colour = "lr") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_12_learning_rate.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Effet de la largeur de la couche cachée
# =============================================================================
widths <- c(1, 2, 3, 5, 10, 20, 50)
tabB <- data.frame(d1 = widths, train = NA, test = NA)
for (i in seq_along(widths)) {
  m <- mlp_fit(tr$x, tr$y, hidden = widths[i], activation = "tanh", loss = "mse",
               epochs = 600, lr = 0.03, batch = 32, seed = 2)
  tabB$train[i] <- mean((predict_mlp(m, tr$x) - tr$y)^2)
  tabB$test[i]  <- mean((predict_mlp(m, te$x) - te$y)^2)
}
cat("\n=== (B) Effet de la largeur cachée d1 (EQM) ===\n")
print(round(tabB, 4), row.names = FALSE)
cat("d1 trop petit : sous-ajustement (biais) ; d1 grand : capacité suffisante.\n")

ggB <- ggplot(reshape(tabB, varying = c("train", "test"), v.names = "eqm",
                      timevar = "ensemble", times = c("train", "test"), direction = "long"),
              aes(d1, eqm, colour = ensemble)) +
  geom_line(linewidth = 1) + geom_point(size = 2) + scale_x_log10() +
  labs(title = "MLP : effet de la largeur de la couche cachée",
       subtitle = "d1 trop petit sous-ajuste ; augmenter d1 réduit le biais",
       x = "largeur cachée d1 (log)", y = "EQM", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_12_largeur.png"), ggB, width = 8, height = 5, dpi = 120)

# =============================================================================
# (C) Précision du gradient numérique vs eps
# =============================================================================
set.seed(5); n <- 30; d0 <- 3; d1 <- 6
Xg <- matrix(rnorm(n * d0), n, d0); Yg <- matrix(rnorm(n), n, 1)
pg <- list(W1 = matrix(rnorm(d0 * d1), d0, d1), b1 = rnorm(d1),
           W2 = matrix(rnorm(d1), d1, 1), b2 = rnorm(1))
ga <- mlp_backward(pg, Xg, Yg, "tanh", "mse")
relerr <- function(a, b) { fa <- unlist(a); fb <- unlist(b)
  sqrt(sum((fa - fb)^2)) / (sqrt(sum(fa^2)) + sqrt(sum(fb^2))) }
epss <- 10^seq(-1, -10, by = -1)
errC <- data.frame(eps = epss,
  err = sapply(epss, function(e) relerr(ga, mlp_numgrad(pg, Xg, Yg, "tanh", "mse", eps = e))))
cat("\n=== (C) Erreur relative gradient numérique vs analytique selon eps ===\n")
print(format(errC, digits = 2, scientific = TRUE), row.names = FALSE)
cat("Décroît en O(eps^2) puis remonte (erreurs d'arrondi) : minimum vers eps~1e-6.\n")

ggC <- ggplot(errC, aes(eps, err)) + geom_line(linewidth = 1, colour = "#2166ac") +
  geom_point(size = 2, colour = "#2166ac") + scale_x_log10() + scale_y_log10() +
  labs(title = "Vérification du gradient : erreur relative selon eps",
       subtitle = "Différences centrées : O(eps^2) puis plancher d'arrondi (minimum ~1e-6)",
       x = expression(epsilon~"(log)"), y = "erreur relative (log)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_12_gradcheck.png"), ggC, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir,
    "/mc_12_learning_rate.png, mc_12_largeur.png, mc_12_gradcheck.png\n", sep = "")
