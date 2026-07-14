# mlfromscratch 0.0.0.9000

Projet de mémoire — réimplémentation en R base de 35 modules (0–34)
d'apprentissage et d'économétrie, validés contre les packages de référence et
étudiés par Monte Carlo.

* Fondations et cœur (Modules 0–16) : algèbre linéaire et optimiseurs avec
  preuves, OLS et inférence, robustesse/GLS, GLM par IRLS, régularisation,
  variables instrumentales, validation, théorie de l'apprentissage, arbres,
  forêts, boosting, M-estimation, non supervisé, MLP, interprétabilité, DML.
* Extensions économétriques (Modules 17–21) : bootstrap, GMM, prédiction
  conforme, régression quantile, panel/effets fixes ; optimiseurs Nesterov et
  L-BFGS.
* Haute dimension et causalité (Modules 22–25) : lasso débiaisé, analyse de
  sensibilité OVB, fonctions d'influence (fil unificateur), DiD à adoption
  échelonnée (Callaway-Sant'Anna, Sun-Abraham).
* Frontières (Modules 26–28) : double descente, processus gaussiens, autodiff.
* Compléments économétriques et statistiques (Modules 29–34) : diagnostics
  (Breusch-Pagan, White, Hausman, Sargan…), variables dépendantes limitées
  (probit, Tobit, Heckman), séries temporelles (AR/ARMA, ADF), non paramétrique
  et RDD, MCMC, méthode delta et tests multiples.
* Contribution originale : diagramme de phase des estimateurs ; chapitre de
  théorie (efficacité semiparamétrique du DML) ; mémoire écrit défendu.
* Infrastructure : compendium reproductible (`renv`, `_targets.R`, `Dockerfile`,
  `CITATION.cff`), site `pkgdown`, bibliographie consolidée, glossaire.
* Tests : plus de 500 vérifications de conformité et de propriétés, exécutées en
  processus isolés (`run_tests.R`), 0 échec.
