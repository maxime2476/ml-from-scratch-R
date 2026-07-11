# =============================================================================
# Monte Carlo — Module 10 : gradient boosting
#  (A) Trajectoires de perte train/test selon M, comparées à gbm (qualitatif).
#  (B) Rôle du taux d'apprentissage nu : petit nu + plus d'arbres généralise
#      mieux ; grand nu sur-ajuste plus vite.
#  (C) Newton vs gradient (log-loss) : le pas de Newton par feuille (éq. 10.7)
#      converge mieux que la feuille = moyenne des pseudo-résidus.
# =============================================================================

for (f in c("08_cart", "10_boosting")) source(file.path("R", paste0(f, ".R")))
suppressMessages({library(ggplot2); library(gbm)})

set.seed(2026)
out_dir <- "simulations/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gen <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$yr <- sin(2 * d$x1) + d$x2 - 0.5 * d$x3 + rnorm(n, sd = 0.5)
  d$yb <- rbinom(n, 1, plogis(1.2 * d$x1 - 0.8 * d$x2 + 0.5 * d$x3))
  d
}
tr <- gen(500, 1); te <- gen(2000, 2)
M <- 300

# =============================================================================
# (A) Trajectoires de perte L2, comparées à gbm
# =============================================================================
fit <- gradient_boost(yr ~ x1 + x2 + x3, tr, "l2", M = M, nu = 0.05,
                      max_depth = 3, min_leaf = 10)
train_path <- boost_loss_path(fit, tr, tr$yr)
test_path  <- boost_loss_path(fit, te, te$yr)
gb <- gbm(yr ~ x1 + x2 + x3, data = tr, distribution = "gaussian", n.trees = M,
          shrinkage = 0.05, interaction.depth = 3, n.minobsinnode = 10, bag.fraction = 1)
gbm_train <- gb$train.error

cat("=== (A) L2 boosting : trajectoires de perte ===\n")
cat(sprintf("train MSE : M=1 -> %.3f ; M=%d -> %.3f\n", train_path[1], M, train_path[M]))
cat(sprintf("test  MSE : M=1 -> %.3f ; M=%d -> %.3f (min %.3f a M=%d)\n",
            test_path[1], M, test_path[M], min(test_path), which.min(test_path)))
cat(sprintf("corrélation trajectoire train (mine vs gbm) = %.4f\n", cor(train_path, gbm_train)))

dfA <- rbind(data.frame(M = 1:M, loss = train_path, courbe = "train (mine)"),
             data.frame(M = 1:M, loss = test_path,  courbe = "test (mine)"),
             data.frame(M = 1:M, loss = gbm_train,  courbe = "train (gbm)"))
ggA <- ggplot(dfA, aes(M, loss, colour = courbe)) + geom_line(linewidth = 0.9) +
  labs(title = "Gradient boosting L2 : trajectoires de perte",
       subtitle = "train décroît sans cesse ; test atteint un minimum puis remonte (sur-ajustement)",
       x = "nombre d'arbres M", y = "EQM", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_10_trajectoire.png"), ggA, width = 8, height = 5, dpi = 120)

# =============================================================================
# (B) Rôle du taux d'apprentissage nu
# =============================================================================
nus <- c(0.01, 0.05, 0.2, 0.5)
dfB <- do.call(rbind, lapply(nus, function(nu) {
  f <- gradient_boost(yr ~ x1 + x2 + x3, tr, "l2", M = M, nu = nu, max_depth = 3, min_leaf = 10)
  data.frame(M = 1:M, test = boost_loss_path(f, te, te$yr), nu = factor(nu))
}))
bestB <- aggregate(test ~ nu, dfB, function(z) c(min = min(z), at = which.min(z)))
cat("\n=== (B) Taux d'apprentissage : test MSE minimal ===\n")
for (i in seq_len(nrow(bestB)))
  cat(sprintf("nu=%-5s : test MSE min = %.3f (a M = %d)\n",
              as.character(bestB$nu[i]), bestB$test[i, "min"], bestB$test[i, "at"]))
cat("Petit nu : min plus bas mais atteint plus tard ; grand nu : sur-ajuste vite.\n")

ggB <- ggplot(dfB, aes(M, test, colour = nu)) + geom_line(linewidth = 0.9) +
  labs(title = "Gradient boosting : effet du taux d'apprentissage",
       subtitle = "Petit nu généralise mieux (min plus bas) mais demande plus d'arbres",
       x = "nombre d'arbres M", y = "EQM de test", colour = "nu") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_10_learning_rate.png"), ggB, width = 8, height = 5, dpi = 120)

# =============================================================================
# (C) Newton vs gradient (log-loss)
# =============================================================================
fn <- gradient_boost(yb ~ x1 + x2 + x3, tr, "logloss", M = M, nu = 0.1,
                     max_depth = 3, min_leaf = 10, newton = TRUE)
fg <- gradient_boost(yb ~ x1 + x2 + x3, tr, "logloss", M = M, nu = 0.1,
                     max_depth = 3, min_leaf = 10, newton = FALSE)
ln <- boost_loss_path(fn, te, te$yb); lg <- boost_loss_path(fg, te, te$yb)
cat("\n=== (C) Log-loss : Newton (éq. 10.7) vs gradient (feuille = moyenne résidus) ===\n")
cat(sprintf("test log-loss a M=%d : Newton=%.4f  gradient=%.4f\n", M, ln[M], lg[M]))
cat(sprintf("test log-loss min    : Newton=%.4f (M=%d)  gradient=%.4f (M=%d)\n",
            min(ln), which.min(ln), min(lg), which.min(lg)))
cat("Le pas de Newton (hessienne) converge plus efficacement que le gradient seul.\n")

dfC <- rbind(data.frame(M = 1:M, loss = ln, methode = "Newton (éq. 10.7)"),
             data.frame(M = 1:M, loss = lg, methode = "gradient (moyenne résidus)"))
ggC <- ggplot(dfC, aes(M, loss, colour = methode)) + geom_line(linewidth = 0.9) +
  labs(title = "Log-loss boosting : Newton vs gradient",
       subtitle = "Le poids de Newton par feuille (pondéré par la hessienne) descend plus vite",
       x = "nombre d'arbres M", y = "log-loss de test", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_10_newton.png"), ggC, width = 8, height = 5, dpi = 120)

cat("\nGraphiques -> ", out_dir,
    "/mc_10_trajectoire.png, mc_10_learning_rate.png, mc_10_newton.png\n", sep = "")
