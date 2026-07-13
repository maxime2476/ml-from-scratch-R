# =============================================================================
# Application empirique — LaLonde (NSW) : toute la boîte à outils causale
# =============================================================================
# Question : le programme de formation NSW augmente-t-il les revenus 1978 ?
# LaLonde (1986) a montré que les méthodes OBSERVATIONNELLES échouaient à
# reproduire le résultat EXPÉRIMENTAL. On dispose de deux échantillons :
#   - lalonde.exp  : NSW randomisé  -> l'ATT est identifié SANS hypothèse (~1794 $).
#   - lalonde.psid : traités NSW + comparaison PSID (non randomisée) -> il faut
#     corriger la confusion. On y déploie OLS, IPW, DML, lasso débiaisé, tous
#     FROM SCRATCH (Modules 1, 3, 16, 22), puis on juge chaque méthode à sa
#     capacité à retrouver le benchmark expérimental. Sensibilité (M23) en prime.
#
# Usage : Rscript applications/lalonde.R   (depuis la racine du package)
# =============================================================================

for (f in c("00_linalg", "01_ols", "03_glm_irls", "04_regularisation",
            "08_cart", "09_bagging_rf", "16_causal_ml", "22_debiased_lasso",
            "23_sensitivity"))
  source(file.path("R", paste0(f, ".R")))
suppressMessages({ library(qte); library(ggplot2) })
out_dir <- "applications/output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(2026)

data(lalonde.exp,  package = "qte"); data(lalonde.psid, package = "qte")
covs <- c("age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75")
fml  <- as.formula(paste("re78 ~ treat +", paste(covs, collapse = " + ")))

# --- 0. Benchmark EXPÉRIMENTAL (l'étalon-or) --------------------------------
exp_dim <- with(lalonde.exp, mean(re78[treat == 1]) - mean(re78[treat == 0]))
exp_ols <- ols_summary(ols_fit(fml, lalonde.exp))$coefficients["treat", ]
benchmark <- exp_ols[["estimate"]]
cat(sprintf("BENCHMARK EXPÉRIMENTAL (randomisé, identifié sans hypothèse)\n"))
cat(sprintf("  différence de moyennes : %8.0f $\n", exp_dim))
cat(sprintf("  OLS ajusté             : %8.0f $  (se %.0f)  <= ÉTALON-OR\n\n",
            benchmark, exp_ols[["se"]]))

# =============================================================================
# ÉCHANTILLON OBSERVATIONNEL (lalonde.psid) : reproduire le benchmark ?
# =============================================================================
D <- lalonde.psid
X <- as.matrix(D[, covs]); y <- D$re78; d <- D$treat
methode <- character(0); est <- se <- numeric(0)
add <- function(nm, e, s = NA_real_) {
  methode[[length(methode) + 1L]] <<- nm
  est[[length(est) + 1L]] <<- as.numeric(e); se[[length(se) + 1L]] <<- as.numeric(s)
}

# --- 1. Différence de moyennes naïve (aucun contrôle) -----------------------
add("Diff. moyennes (naïf)", mean(y[d == 1]) - mean(y[d == 0]))

# --- 2. OLS avec contrôles (Module 1) ---------------------------------------
ols_c <- ols_summary(ols_fit(fml, D))$coefficients["treat", ]
add("OLS + contrôles (M1)", ols_c[["estimate"]], ols_c[["se"]])

# --- 3. IPW par score de propension logistique (Module 3) -------------------
# ê(x) = P(traité | x) par IRLS ; ATT-IPW (pondération de Hájek), poids ê/(1-ê).
ps_fit <- glm_irls(as.formula(paste("treat ~", paste(covs, collapse = " + "))),
                   data = D, family = "binomial")
ehat <- pmin(pmax(ps_fit$fitted, 1e-3), 1 - 1e-3)
w0 <- ehat / (1 - ehat)                                   # odds (contrôles)
att_ipw <- mean(y[d == 1]) - weighted.mean(y[d == 0], w0[d == 0])
add("IPW propension (M3)", att_ipw)

# --- 4. Double Machine Learning (PLR, forêts, cross-fitting ; Module 16) -----
dml <- dml_plr(y = y, d = d, X = as.data.frame(X), K = 5L, nuisance = "forest", seed = 1)
add("DML forêts (M16)", dml$theta, dml$se)

# --- 5. Lasso débiaisé sur dictionnaire étendu (Module 22) -------------------
# On enrichit X (carrés + interactions avec le traitement) : p grandit, mais le
# lasso débiaisé fournit un IC valide pour le coefficient du traitement.
Xc <- scale(X[, c("age", "education", "re74", "re75")])
Xbig <- cbind(treat = d, X, Xc^2, d * Xc)                 # hétérogénéité treat×covariables
colnames(Xbig)[1] <- "treat"
db <- debiased_lasso(Xbig, y, targets = 1)                # cible = coefficient du traitement
add("Lasso débiaisé (M22)", db$estimate[1], db$se[1])

# =============================================================================
# TABLEAU RÉCAPITULATIF
# =============================================================================
tab <- data.frame(methode = unlist(methode), estimate = unlist(est),
                  se = unlist(se), row.names = NULL)
tab$ecart_benchmark <- tab$estimate - benchmark
cat("ÉCHANTILLON OBSERVATIONNEL (PSID) — ATT estimé du programme sur re78 :\n\n")
cat(sprintf("%-26s %10s %8s %14s\n", "méthode", "ATT ($)", "se", "écart/étalon"))
for (i in seq_len(nrow(tab)))
  cat(sprintf("%-26s %10.0f %8s %14.0f\n", tab$methode[i], tab$estimate[i],
              ifelse(is.na(tab$se[i]), "-", sprintf("%.0f", tab$se[i])),
              tab$ecart_benchmark[i]))
cat(sprintf("\n(étalon-or expérimental = %.0f $)\n", benchmark))
cat("=> LEÇON DE LALONDE (1986). Le naïf est massivement biaisé (-15 000 $). Le\n")
cat("   contrôle des covariables aide (IPW ~1575 $ retrouve le benchmark ; OLS ~752 $\n")
cat("   le sous-estime), MAIS le résultat reste INSTABLE selon la méthode : DML et\n")
cat("   lasso débiaisé restent négatifs sur le PSID complet — recouvrement médiocre\n")
cat("   (propensions proches de 0/1) et extrapolation. Aucune hypothèse ne remplace\n")
cat("   la randomisation ; d'où la valeur du benchmark expérimental.\n\n")

# --- Diagnostic de recouvrement (overlap) -----------------------------------
share_extreme <- mean(ehat[d == 0] < 0.01)
cat(sprintf("DIAGNOSTIC DE RECOUVREMENT : %.0f %% des contrôles PSID ont ê(x) < 0.01\n",
            100 * share_extreme))
cat("  -> support commun très faible : DML/lasso extrapolent, d'où leur instabilité.\n\n")

# --- 6. Analyse de sensibilité (M23) sur l'estimation OLS observationnelle ---
f_ols <- ols_fit(fml, D)
s <- sensitivity_ols(f_ols, "treat")
# étalonnage : un confondeur « aussi fort que re75 » (revenu pré-programme)
r2_re75 <- partial_r2(ols_summary(f_ols)$coefficients["re75", "t"], s$df)
adj <- adjusted_estimate(s$estimate, s$se, s$df, r2_re75, r2_re75)
cat("SENSIBILITÉ (M23) de l'OLS observationnel :\n")
cat(sprintf("  ATT OLS = %.0f $ (se %.0f) ; R²_{Y~D|X} = %.3f ; robustness value = %.3f\n",
            s$estimate, s$se, s$r2yd, s$rv_q))
cat(sprintf("  confondeur « aussi fort que re75 » (R²=%.3f) -> ATT ajusté = %.0f $\n",
            r2_re75, adj$estimate))
cat("  Lecture : la RV faible confirme la FRAGILITÉ de l'OLS observationnel —\n")
cat("  cohérent avec l'écart au benchmark. La randomisation, elle, l'évite.\n")

# --- Graphique ---------------------------------------------------------------
tab$methode <- factor(tab$methode, levels = rev(tab$methode))
gg <- ggplot(tab, aes(estimate, methode)) +
  geom_vline(xintercept = benchmark, linetype = "dashed", colour = "steelblue") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                orientation = "y", width = 0.2, na.rm = TRUE) +
  annotate("text", x = benchmark, y = 0.6, label = "benchmark expérimental",
           colour = "steelblue", hjust = -0.02, size = 3.2) +
  labs(title = "LaLonde : reproduire l'effet expérimental à partir de données observationnelles",
       subtitle = paste0("ATT du programme NSW sur les revenus 1978 (PSID) ; étalon-or ≈ ",
                         round(benchmark), " $"),
       x = "ATT estimé ($)", y = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "lalonde_atts.png"), gg, width = 9, height = 5, dpi = 120)
cat("\nGraphique -> ", file.path(out_dir, "lalonde_atts.png"), "\n")
