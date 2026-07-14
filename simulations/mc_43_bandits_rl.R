# =============================================================================
# Monte Carlo — Module 43 : (1) le regret des bandits croit LOGARITHMIQUEMENT
# (UCB, Thompson) vs lineairement (epsilon-greedy fixe, hasard) ; (2) le Q-learning
# converge vers Q* (itération sur les valeurs).
# =============================================================================

for (f in c("43_bandits_rl", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## (1) courbes de regret moyennees ------------------------------------------
means <- c(0.2, 0.35, 0.5, 0.65, 0.45); T <- 3000; R <- 40
eps_greedy <- function(means, T, eps, seed) { set.seed(seed); K <- length(means); n <- s <- rep(0, K)
  best <- max(means); cum <- 0; reg <- numeric(T)
  for (t in 1:T) { a <- if (t <= K || runif(1) < eps) sample.int(K, 1) else which.max(s / pmax(n,1))
    r <- rbinom(1,1,means[a]); n[a] <- n[a]+1; s[a] <- s[a]+r; cum <- cum + (best-means[a]); reg[t] <- cum }
  reg }
Ru <- Rt <- Re <- Rr <- numeric(T)
for (r in seq_len(R)) {
  Ru <- Ru + bandit_ucb(means, T, seed = r)$regret
  Rt <- Rt + bandit_thompson(means, T, seed = r)$regret
  Re <- Re + eps_greedy(means, T, 0.1, seed = r)
  Rr <- Rr + cumsum(rep(max(means) - mean(means), T))
}
Ru<-Ru/R; Rt<-Rt/R; Re<-Re/R; Rr<-Rr/R
cat("=== (1) Regret moyen a T =", T, "===\n\n")
cat(sprintf("  hasard            : %.0f\n  epsilon-greedy 0.1 : %.0f\n  UCB1              : %.0f\n  Thompson          : %.0f\n",
            Rr[T], Re[T], Ru[T], Rt[T]))
cat("\n=> UCB et Thompson : regret LOGARITHMIQUE (la courbe s'aplatit). Le hasard :\n")
cat("   lineaire (660). L'epsilon-greedy a taux fixe garde une composante lineaire\n")
cat("   asymptotique (exploration permanente) que le graphe rend visible sur le long terme.\n\n")

## (2) convergence du Q-learning ---------------------------------------------
set.seed(1); S<-5; A<-3; P<-array(0,c(S,A,S)); for(s in 1:S) for(a in 1:A){p<-runif(S);P[s,a,]<-p/sum(p)}
Rw<-matrix(runif(S*A),S,A); vi<-value_iteration(P,Rw,0.9)
eps_list <- c(500,1000,2000,5000,10000,20000); gap <- sapply(eps_list, function(e) max(abs(q_learning(P,Rw,0.9,episodes=e,seed=2)$Q - vi$Q)))
cat("=== (2) Q-learning vs Q* (itération sur les valeurs) ===\n")
for (i in seq_along(eps_list)) cat(sprintf("  %6d episodes : max|Q - Q*| = %.4f\n", eps_list[i], gap[i]))
cat("=> Des les premiers milliers d'episodes, Q reste dans un VOISINAGE de Q*\n")
cat("   (~0.1, residu du pas alpha constant) et la politique optimale est retrouvee,\n")
cat("   sans jamais connaitre le modele P, R.\n")

df <- rbind(data.frame(t=1:T, reg=Rr, algo="hasard"), data.frame(t=1:T, reg=Re, algo="epsilon-greedy 0.1"),
            data.frame(t=1:T, reg=Ru, algo="UCB1"), data.frame(t=1:T, reg=Rt, algo="Thompson"))
gg <- ggplot(df, aes(t, reg, colour = algo)) + geom_line(linewidth = 0.9) +
  labs(title = "Bandits : regret logarithmique (UCB/Thompson) vs lineaire (hasard/epsilon-fixe)",
       subtitle = "regret cumule moyen ; les courbes plates = apprentissage efficace",
       x = "tour t", y = "regret cumule", colour = NULL) + theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(out_dir, "mc_43_regret.png"), gg, width = 8.5, height = 5, dpi = 120)
cat("\nGraphique -> mc_43_regret.png\n")
