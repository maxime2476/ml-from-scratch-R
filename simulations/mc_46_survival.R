# =============================================================================
# Monte Carlo — Module 46 : (1) le Cox recupere le vrai hazard ratio (taille/
# puissance du test), (2) ignorer la CENSURE (OLS sur les temps observes) biaise
# gravement, la ou Kaplan-Meier/Cox la gerent correctement.
# =============================================================================

for (f in c("46_survival", "mc_tools")) source(file.path("R", paste0(f, ".R")))
suppressMessages(library(ggplot2))
set.seed(2026)
out_dir <- "simulations/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gen <- function(n, beta) { x <- rnorm(n); tt <- rexp(n, rate = exp(beta * x) * 0.1)
  cens <- runif(n, 0, 12); list(time = pmin(tt, cens), event = as.integer(tt <= cens), x = x, tt = tt) }

## (1) Cox : recuperation et test du hazard ratio ----------------------------
R <- 300; n <- 300; z <- qnorm(0.975); bias <- rej0 <- rej1 <- 0
for (r in seq_len(R)) {
  d0 <- gen(n, 0); d1 <- gen(n, 0.8)
  c1 <- cox_ph(d1$time, d1$event, d1$x); c0 <- cox_ph(d0$time, d0$event, d0$x)
  bias <- bias + (as.numeric(c1$coefficients) - 0.8)
  rej0 <- rej0 + (abs(c0$coefficients / c0$se) > z)          # H0 : beta=0
  rej1 <- rej1 + (abs(c1$coefficients / c1$se) > z)          # H1 : beta=0.8
}
cat(sprintf("=== (1) Cox (vrai beta = 0.8, n=%d, R=%d) ===\n\n", n, R))
cat(sprintf("  biais moyen de beta_hat : %+.3f\n", bias / R))
cat(sprintf("  taille du test (beta=0) : %.3f  | puissance (beta=0.8) : %.3f\n\n", rej0 / R, rej1 / R))

## (2) Ignorer la censure biaise l'estimation de la survie mediane -----------
set.seed(1); d <- gen(2000, 0)
km <- kaplan_meier(d$time, d$event)
med_km <- km$time[which.min(abs(km$surv - 0.5))]            # mediane de survie (KM)
med_naif <- median(d$time)                                  # mediane des temps OBSERVES (censure incluse)
med_vrai <- median(d$tt)                                    # vraie mediane (temps complets)
cat("=== (2) Mediane de survie (censure ~ ", round(mean(d$event == 0) * 100), "% ) ===\n", sep = "")
cat(sprintf("  vraie                          : %.2f\n", med_vrai))
cat(sprintf("  Kaplan-Meier (gere la censure) : %.2f\n", med_km))
cat(sprintf("  naif (mediane des observes)    : %.2f  (sous-estime : la censure tronque)\n", med_naif))
cat("\n=> Le Cox recupere le hazard ratio ; Kaplan-Meier corrige le biais de censure\n")
cat("   que l'analyse naive des durees observees introduit.\n")

# figure : courbes KM par groupe
d2 <- gen(600, 1.0); grp <- factor(ifelse(d2$x > 0, "risque eleve", "risque faible"))
df <- do.call(rbind, lapply(levels(grp), function(g) {
  k <- kaplan_meier(d2$time[grp == g], d2$event[grp == g]); data.frame(t = c(0, k$time), S = c(1, k$surv), groupe = g) }))
gg <- ggplot(df, aes(t, S, colour = groupe)) + geom_step(linewidth = 1) +
  labs(title = "Courbes de survie de Kaplan-Meier par groupe de risque",
       subtitle = "le groupe a covariable elevee (hazard ratio > 1) survit moins longtemps",
       x = "temps", y = "survie S(t)", colour = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "mc_46_survival.png"), gg, width = 8, height = 5, dpi = 120)
cat("\nGraphique -> mc_46_survival.png\n")
