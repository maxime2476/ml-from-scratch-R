# Changelog

## mlfromscratch 0.0.0.9000

Projet de mémoire — réimplémentation en R base de 48 modules (0–47)
d’apprentissage et d’économétrie, validés contre les packages de
référence et étudiés par Monte Carlo.

- Fondations et cœur (Modules 0–16) : algèbre linéaire et optimiseurs
  avec preuves, OLS et inférence, robustesse/GLS, GLM par IRLS,
  régularisation, variables instrumentales, validation, théorie de
  l’apprentissage, arbres, forêts, boosting, M-estimation, non
  supervisé, MLP, interprétabilité, DML.
- Extensions économétriques (Modules 17–21) : bootstrap, GMM, prédiction
  conforme, régression quantile, panel/effets fixes ; optimiseurs
  Nesterov et L-BFGS.
- Haute dimension et causalité (Modules 22–25) : lasso débiaisé, analyse
  de sensibilité OVB, fonctions d’influence (fil unificateur), DiD à
  adoption échelonnée (Callaway-Sant’Anna, Sun-Abraham).
- Frontières (Modules 26–28) : double descente, processus gaussiens,
  autodiff.
- Compléments économétriques et statistiques (Modules 29–34) :
  diagnostics (Breusch-Pagan, White, Hausman, Sargan…), variables
  dépendantes limitées (probit, Tobit, Heckman), séries temporelles
  (AR/ARMA, ADF), non paramétrique et RDD, MCMC, méthode delta et tests
  multiples.
- Réseaux de neurones profonds (Modules 35–38) : optimiseurs (Adam) et
  régularisation (dropout, batch norm), CNN, RNN/LSTM, attention et
  Transformer.
- Machine learning — compléments (Modules 39–43) : SVM (dual QP, noyau),
  classifieurs génératifs (LDA/QDA/Naive Bayes), clustering avancé
  (hiérarchique, DBSCAN, spectral), réduction de dimension (kernel PCA,
  ICA, NMF, t-SNE), bandits et apprentissage par renforcement (UCB,
  Thompson, Q-learning).
- Économétrie avancée (Modules 44–47) : VAR (Granger, IRF), GARCH,
  analyse de survie (Kaplan-Meier, Cox), panel avancé (contrôle
  synthétique, panel dynamique).
- Contribution originale : diagramme de phase des estimateurs ; chapitre
  de théorie (efficacité semiparamétrique du DML) ; mémoire écrit
  défendu.
- Infrastructure : compendium reproductible (`renv`, `_targets.R`,
  `Dockerfile`, `CITATION.cff`), site `pkgdown`, bibliographie
  consolidée, glossaire.
- Tests : 622 vérifications de conformité et de propriétés, exécutées en
  processus isolés (`run_tests.R`), 0 échec.
