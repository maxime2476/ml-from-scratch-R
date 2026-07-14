# =============================================================================
# Monte Carlo — Module 25 : le TWFE échoue, Callaway-Sant'Anna et Sun-Abraham
# récupèrent l'ATT sous adoption échelonnée à effets dynamiques hétérogènes.
# =============================================================================

source(file.path("R", "25_did.R"))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N <- 400; Tt <- 8; cohorts <- c(3, 5, 7, Inf)
gen <- function(seed) {
  set.seed(seed)
  id <- rep(1:N, each = Tt); t <- rep(1:Tt, N)
  g <- rep(sample(cohorts, N, replace = TRUE), each = Tt)
  ai <- rep(rnorm(N), each = Tt); gt <- rep(0.3 * (1:Tt), N)
  D <- as.integer(t >= g); expo <- pmax(t - g + 1, 0); expo[is.infinite(g)] <- 0
  tau <- ifelse(D == 1, 1 + 0.5 * expo + ifelse(g == 3, 1, 0), 0)   # dynamique + hétéro
  list(data = data.frame(id, t, g, y = ai + gt + tau + rnorm(length(id))),
       att_true = mean(tau[D == 1]))
}

R <- 400; b_tw <- b_cs <- b_sa <- att_true <- numeric(R)
for (r in seq_len(R)) {
  s <- gen(r); d <- s$data; att_true[r] <- s$att_true
  b_tw[r] <- twfe(d)$coef
  b_cs[r] <- aggregate_att(att_gt(d, control = "never"), "simple")
  b_sa[r] <- sunab(d)$att
}
cat("=== ATT sous adoption échelonnée (effets dynamiques hétérogènes) ===\n\n")
cat(sprintf("  ATT vrai (moyenne)          : %.3f\n", mean(att_true)))
cat(sprintf("  TWFE            : %.3f  (biais %+.3f)\n", mean(b_tw), mean(b_tw - att_true)))
cat(sprintf("  Callaway-Sant'Anna : %.3f  (biais %+.3f)\n", mean(b_cs), mean(b_cs - att_true)))
cat(sprintf("  Sun-Abraham     : %.3f  (biais %+.3f)\n", mean(b_sa), mean(b_sa - att_true)))
cat("\n=> Le TWFE sous-estime gravement (comparaisons interdites, poids négatifs) ;\n")
cat("   Callaway-Sant'Anna et Sun-Abraham récupèrent l'ATT vrai.\n")

# --- Figure 1 : distribution des estimateurs vs vérité ----------------------
df <- data.frame(est = c(b_tw, b_cs, b_sa),
                 methode = rep(c("TWFE", "Callaway-Sant'Anna", "Sun-Abraham"), each = R))
df$methode <- factor(df$methode, levels = c("TWFE", "Callaway-Sant'Anna", "Sun-Abraham"))
g1 <- ggplot(df, aes(est, fill = methode)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = mean(att_true), linetype = "dashed") +
  annotate("text", x = mean(att_true), y = 0, label = "ATT vrai", vjust = -0.5, hjust = -0.05) +
  labs(title = "DiD échelonné : le TWFE est biaisé, CS/SA récupèrent l'ATT",
       subtitle = paste0("distribution des estimateurs sur ", R, " simulations"),
       x = "ATT estimé", y = "densité", fill = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_25_did.png"), g1, width = 8, height = 5, dpi = 120)

# --- Figure 2 : event-study Sun-Abraham vs vérité (un tirage) ----------------
s <- gen(1); d <- s$data
sa <- sunab(d)$es
truth_es <- data.frame(e = 0:4, att = 1 + 0.5 * (1:5))     # profil dynamique moyen (approx.)
es <- merge(sa, transform(truth_es, vrai = att)[, c("e", "vrai")], by = "e", all.x = TRUE)
g2 <- ggplot(es[es$e >= -3, ], aes(e, att)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = -0.5, linetype = "dotted") +
  geom_line(colour = "#00798c") + geom_point(colour = "#00798c") +
  geom_point(aes(y = vrai), colour = "#d1495b", shape = 4, size = 2, na.rm = TRUE) +
  labs(title = "Event-study Sun-Abraham (points) vs effet dynamique vrai (croix)",
       subtitle = "les pré-périodes (e < 0) sont à ~0 : pas de fausse tendance",
       x = "ancienneté du traitement (e = t - g)", y = "effet") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_25_eventstudy.png"), g2, width = 8, height = 5, dpi = 120)
cat("\nGraphiques -> mc_25_did.png, mc_25_eventstudy.png\n")
