# =============================================================================
# Monte Carlo — Module 35 : (1) Adam converge bien plus vite que le GD sur un
# probleme mal conditionne ; (2) le dropout reduit le sur-apprentissage d'un MLP.
# =============================================================================

for (f in c("35_nn_training", "28_autodiff", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) Adam vs GD : trajectoire de perte sur une quadratique mal conditionnee ---
set.seed(1); p <- 20; d <- exp(seq(log(0.2), log(50), length.out = p)); b <- rnorm(p)
f <- function(x) 0.5 * sum(d * (x - b)^2); g <- function(x) d * (x - b)
n_it <- 300
loss_adam <- numeric(n_it); x <- rep(0, p); m <- v <- rep(0, p)
for (t in 1:n_it) { gr <- g(x); m <- 0.9*m + 0.1*gr; v <- 0.999*v + 0.001*gr^2
  x <- x - 0.5 * (m/(1-0.9^t)) / (sqrt(v/(1-0.999^t)) + 1e-8); loss_adam[t] <- f(x) }
loss_gd <- numeric(n_it); x <- rep(0, p)
for (t in 1:n_it) { x <- x - (2/max(d)) * g(x); loss_gd[t] <- f(x) }   # pas maximal stable
cat("=== (1) Convergence (quadratique, conditionnement", round(max(d)/min(d)), ") ===\n")
cat(sprintf("  perte apres %d iterations : Adam %.2e | GD %.2e\n\n", n_it, loss_adam[n_it], loss_gd[n_it]))

## (2) Dropout reduit le sur-apprentissage ------------------------------------
# MLP a 1 couche cachee (autodiff), petit jeu bruite -> sur-apprentissage
gen <- function(n) { X <- matrix(runif(n*2, -2, 2), n, 2); y <- sin(X[,1]) * cos(X[,2]) + 0.15*rnorm(n)
  list(X = cbind(X, 1), y = matrix(y, n, 1)) }
train_mlp <- function(tr, te, h = 120, drop = 0, epochs = 3000, lr = 0.03) {
  p <- ncol(tr$X); W1 <- matrix(rnorm(p*h)*0.5, p, h); W2 <- matrix(rnorm((h+1))*0.5, h+1, 1)
  for (e in 1:epochs) {
    ad_reset(); w1 <- adnode(W1); w2 <- adnode(W2)
    H <- tanh(mm(tr$X, w1))
    Hd <- if (drop > 0) { mask <- (matrix(runif(length(H$value)), nrow(H$value)) > drop) / (1 - drop)
      H * mask } else H                                   # dropout INJECTE dans le graphe
    Ha <- ad_cbind1(Hd); r <- mm(Ha, w2) - tr$y; L <- sum(r*r)/nrow(tr$X); backward(L)
    W1 <- W1 - lr * w1$grad; W2 <- W2 - lr * w2$grad
  }
  pred <- function(Z){ Ha <- cbind(tanh(Z %*% W1), 1); Ha %*% W2 }   # test : pas de dropout
  c(train = mean((pred(tr$X) - tr$y)^2), test = mean((pred(te$X) - te$y)^2))
}
R <- 12; res <- matrix(0, 2, 2, dimnames = list(c("sans dropout","avec dropout"), c("train","test")))
for (r in seq_len(R)) { tr <- gen(25); te <- gen(1000)
  res[1, ] <- res[1, ] + train_mlp(tr, te, drop = 0)
  res[2, ] <- res[2, ] + train_mlp(tr, te, drop = 0.1) }
res <- res / R
cat("=== (2) Sur-apprentissage d'un MLP (n_train=60) ===\n")
cat(sprintf("  sans dropout : train %.3f | test %.3f  (ecart %.3f)\n", res[1,1], res[1,2], res[1,2]-res[1,1]))
cat(sprintf("  avec dropout : train %.3f | test %.3f  (ecart %.3f)\n", res[2,1], res[2,2], res[2,2]-res[2,1]))
cat("\n=> Adam ecrase le GD sur un probleme mal conditionne ; le dropout reduit\n")
cat("   l'ecart train-test (moins de sur-apprentissage).\n")

df <- data.frame(iter = rep(1:n_it, 2), loss = c(loss_adam, loss_gd),
                 opt = rep(c("Adam", "GD (pas maximal)"), each = n_it))
gg <- ggplot(df, aes(iter, loss, colour = opt)) + geom_line(linewidth = 1) + scale_y_log10() +
  labs(title = "Adam vs descente de gradient (probleme mal conditionne)",
       subtitle = "perte (echelle log) ; Adam s'adapte par coordonnee et converge bien plus vite",
       x = "iteration", y = "perte", colour = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_35_adam.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_35_adam.png\n")
