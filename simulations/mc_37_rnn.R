# =============================================================================
# Monte Carlo — Module 37 : le gradient qui s'evanouit (RNN) vs l'autoroute (LSTM)
# On mesure la SENSIBILITE de l'etat final au PREMIER pas d'entree, en fonction de
# la longueur T de la sequence. Pour le RNN simple (tanh), elle decroit
# exponentiellement (memoire courte) ; pour la LSTM avec porte d'oubli ouverte,
# elle persiste (l'etat de cellule preserve le signal lointain).
# =============================================================================

for (f in c("37_rnn", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages({ library(ggplot2); library(numDeriv) })
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- 3; H <- 5; Ts <- c(5, 10, 20, 40, 60)
sens_rnn <- sens_lstm <- numeric(length(Ts))

for (k in seq_along(Ts)) {
  Tn <- Ts[k]; set.seed(1)
  X <- matrix(rnorm(Tn * d), Tn, d)
  Wxh <- matrix(rnorm(H * d) * 0.5, H, d); Whh <- matrix(rnorm(H * H) * 0.4, H, H)   # rayon spectral < 1
  Why <- diag(H)[1, , drop = FALSE]; bh <- rep(0, H)
  # sensibilite de h_T au premier vecteur d'entree x_1 (norme du jacobien)
  hT_rnn <- function(x1) { Xm <- X; Xm[1, ] <- x1
    rnn_forward(Xm, Wxh, Whh, Why, bh, 0)$H[Tn, ] }
  J_rnn <- jacobian(hT_rnn, X[1, ]); sens_rnn[k] <- sqrt(sum(J_rnn^2))

  W <- matrix(rnorm(H * (d + H)) * 0.4, H, d + H)
  hT_lstm <- function(x1) { Xm <- X; Xm[1, ] <- x1
    lstm_forward(Xm, W, W, W, W, rep(0, H), rep(3, H), rep(0, H), rep(0, H))$H[Tn, ] }  # bf=3 : oubli ouvert
  J_lstm <- jacobian(hT_lstm, X[1, ]); sens_lstm[k] <- sqrt(sum(J_lstm^2))
}

cat("=== Sensibilite de l'etat final au premier pas, selon T ===\n\n")
cat(sprintf("%5s %14s %14s\n", "T", "RNN", "LSTM"))
for (k in seq_along(Ts)) cat(sprintf("%5d %14.2e %14.2e\n", Ts[k], sens_rnn[k], sens_lstm[k]))
cat(sprintf("\n=> RNN : chute d'un facteur %.0e entre T=%d et T=%d (gradient evanoui).\n",
            sens_rnn[1] / sens_rnn[length(Ts)], Ts[1], max(Ts)))
cat("   LSTM : la sensibilite au passe lointain PERSISTE (autoroute de memoire).\n")

df <- rbind(data.frame(T = Ts, sens = sens_rnn, modele = "RNN simple"),
            data.frame(T = Ts, sens = sens_lstm, modele = "LSTM"))
gg <- ggplot(df, aes(T, sens, colour = modele)) + geom_line(linewidth = 1) + geom_point() +
  scale_y_log10() +
  labs(title = "Gradient dans le temps : le RNN oublie, la LSTM se souvient",
       subtitle = "sensibilite (norme du jacobien) de l'etat final au premier pas, echelle log",
       x = "longueur de sequence T", y = "sensibilite au passe lointain", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_37_vanishing.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_37_vanishing.png\n")
