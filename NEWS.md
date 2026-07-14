# mlfromscratch (development version)

Projet de mémoire — réimplémentation en R base de 29 modules d'apprentissage et
d'économétrie, validés contre les packages de référence et étudiés par Monte
Carlo. Historique des grandes étapes.

## Fondations et cœur (Modules 0–16)

- Algèbre linéaire et optimiseurs avec preuves de convergence (M0) ; OLS et
  inférence (M1) ; robustesse/GLS (M2) ; GLM par IRLS (M3) ; régularisation
  ridge/lasso (M4) ; variables instrumentales/2SLS (M5) ; validation (M6).
- Théorie de l'apprentissage (M13) ; modèles non linéaires : KNN, CART, forêts,
  boosting (M7–M10) ; M-estimation (M14) ; non supervisé et MLP (M11–M12) ;
  interprétabilité (M15) ; pont ML ↔ causalité, DML (M16).

## Extensions économétriques (Modules 17–21)

- Bootstrap (M17) ; GMM (M18) ; prédiction conforme (M19) ; régression quantile
  (M20) ; panel / effets fixes (M21). Optimiseurs Nesterov et L-BFGS ajoutés (M0).

## Inférence en haute dimension et causalité (Modules 22–25)

- Lasso débiaisé (M22) ; analyse de sensibilité OVB (M23) ; **fonctions
  d'influence**, le fil unificateur (M24) ; **DiD à adoption échelonnée** —
  poids négatifs, Callaway-Sant'Anna, Sun-Abraham (M25).

## Frontières (Modules 26–28)

- Double descente (M26) ; processus gaussiens et noyaux (M27) ; différentiation
  automatique en mode inverse (M28).

## Recherche et infrastructure

- Contribution originale : **diagramme de phase des estimateurs** (carte de
  décision « quelle méthode quand »).
- Chapitre de théorie approfondie : preuve de l'efficacité semiparamétrique du DML.
- Études : LaLonde, double sélection, benchmark from-scratch vs référence.
- Application Shiny interactive ; mémoire écrit défendu et livre Quarto.
- Compendium reproductible : `renv`, `_targets.R`, `Dockerfile`, `CITATION.cff`,
  site `pkgdown`.
- Tests : **375+ vérifications** de conformité et de **propriétés**, exécutées en
  processus isolés (`run_tests.R`), 0 échec.
