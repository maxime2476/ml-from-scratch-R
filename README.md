# ML from Scratch in R

**Déconstruction mathématique et statistique des modèles d'apprentissage automatique**

Projet de fin d'études (statistique / économétrie). Chaque modèle de *machine
learning* est réimplémenté en **R base**, à partir de sa dérivation
mathématique complète, puis validé numériquement contre les packages de
référence et étudié par simulation Monte Carlo.

## Objectif

L'objectif n'est **pas** de produire du code performant, mais de **démontrer une
compréhension profonde des fondements** mathématiques, statistiques et
économétriques des modèles. Chaque implémentation est adossée à :

1. une **dérivation** rédigée depuis les premiers principes (vraisemblance,
   moindres carrés, conditions d'optimalité) ;
2. un **code R base** lisible, dont chaque fonction renvoie explicitement à
   l'équation qu'elle implémente ;
3. une **validation numérique** contre le package de référence (tolérance
   cible `1e-8`, ou justification documentée quand elle n'est pas atteignable) ;
4. une **étude Monte Carlo** avec DGP connu, mesurant biais, variance,
   couverture des intervalles de confiance et puissance des tests.

## Principes méthodologiques

- **Les mathématiques d'abord.** Aucune ligne de code sans dérivation
  préalable.
- **R base pour le cœur des modèles.** Les packages de ML (`glmnet`, `rpart`,
  `caret`, `tidymodels`, …) servent **uniquement** à la validation, jamais à
  l'implémentation.
- **Validation systématique** contre les fonctions de référence.
- **Approche par simulation** : DGP explicite, paramètres vrais fixés,
  ≥ 1000 réplications Monte Carlo.
- **Rigueur inférentielle** : tout estimateur est présenté avec ses propriétés
  (biais, variance, loi asymptotique) et les hypothèses qui les fondent.
- **Auto-suffisance par module.** Chaque module réimplémente localement les
  briques dont il a besoin (notamment son optimiseur), plutôt que d'importer
  celles du Module 0. Le lecteur d'un module n'a jamais à remonter ailleurs pour
  comprendre son code. Le Module 0 reste le lieu **unique** des preuves de
  convergence et des versions canoniques de référence.

## Structure

```
ml-from-scratch-R/
├── README.md                 présentation, table des matières
├── DESCRIPTION               métadonnées du package R minimal
├── R/                        implémentations (R base)
├── derivations/              une dérivation Quarto (.qmd) par module
├── tests/testthat/           tests de conformité aux packages de référence
├── simulations/              études Monte Carlo (un script par question)
└── rapport/                  rapport final Quarto assemblant l'ensemble
```

Pour chaque module `XX`, quatre livrables produits **dans cet ordre** :

| Livrable | Emplacement | Contenu |
|----------|-------------|---------|
| (a) Dérivation | `derivations/XX.qmd` | hypothèses numérotées, dérivation pas à pas, propriétés, encadré « de la math au code » |
| (b) Implémentation | `R/XX.R` | fonctions R base, docstrings roxygen2 renvoyant aux équations |
| (c) Tests | `tests/testthat/test-XX.R` | conformité aux packages de référence (`tolerance = 1e-8`) |
| (d) Monte Carlo | `simulations/mc_XX.R` | DGP, ≥ 1000 réplications, tableaux + graphiques |

## Table des matières (modules)

Les modules sont traités dans l'**ordre de la colonne « # »** (pédagogique),
qui diffère de la numérotation des fichiers : les modules 13 à 16, ajoutés
après coup, s'insèrent à leur place logique. Ordre global :

> 0 → 1 → 2 → 3 → 4 → 5 → 6 → **13** → 7 → 8 → 9 → 10 → **14** → 11 → 12 → **15** → **16**

| Ordre | Fichier | Module | Points clés de dérivation | Référence de validation |
|:-----:|:-------:|--------|---------------------------|-------------------------|
| 1 | `00` | Algèbre linéaire & optimiseurs | QR (Householder), Cholesky, **SVD** (rang, pseudo-inverse Moore-Penrose) ; κ(XᵀX)=κ(X)² ; **catalogue de référence des optimiseurs** : descente de gradient (conv. O(1/k)), Newton (conv. quadratique), coordinate descent, SGD | `qr`, `chol`, `svd` |
| 2 | `01` | OLS et inférence | Gauss-Markov (preuve), lois t/F, Frisch-Waugh-Lovell | `lm`, `summary.lm` |
| 3 | `02` | Hétéroscédasticité & robustesse | sandwich HC0–HC3, Newey-West, WLS/GLS | `sandwich`, `nlme::gls` |
| 4 | `03` | GLM et IRLS | famille exponentielle, Newton ⇔ IRLS, Wald/LR/score | `glm` |
| 5 | `04` | Régularisation | ridge (biais/variance via SVD), lasso (soft-thresholding, coordinate descent réimplémenté localement) | `MASS::lm.ridge`, `glmnet` |
| 6 | `05` | Variables instrumentales | biais d'endogénéité, 2SLS, instruments faibles | `AER::ivreg` |
| 7 | `06` | Validation de modèles | biais-variance (preuve), LOOCV via hat matrix, AIC/BIC | — |
| 8 | `13` | **Théorie de l'apprentissage statistique** | PAC/MRE, Hoeffding (preuve), VC-dim & Sauer-Shelah, complexité de Rademacher, lien avec la régularisation | *illustrations numériques* |
| 9 | `07` | KNN & fléau de la dimension | estimateur local, concentration des distances | `class::knn` |
| 10 | `08` | CART | Gini/entropie, cost-complexity pruning | `rpart` |
| 11 | `09` | Bagging & forêts aléatoires | variance d'estimateurs corrélés, erreur OOB | `randomForest` |
| 12 | `10` | Gradient boosting | descente de gradient fonctionnelle, Newton boosting | `gbm` |
| 13 | `14` | **Le regard économétrique sur le ML** | M-estimation unificatrice (sandwich A⁻¹BA⁻¹), lecture bayésienne/MAP, inférence post-sélection cassée | *illustrations numériques* |
| 14 | `11` | Non supervisé | PCA (2 voies), k-means, EM gaussien | `prcomp`, `mclust` |
| 15 | `12` | MLP minimal | backpropagation, vérification du gradient (SGD réimplémenté localement) | — |
| 16 | `15` | **Interprétabilité post-hoc** | PDP/ICE, valeurs de Shapley (axiomes + preuve), SHAP, forme fermée linéaire | `iml` / `fastshap` |
| 17 | `16` | **Pont ML ↔ causalité** *(module final)* | résultats potentiels, score orthogonal de Neyman, DML + cross-fitting, forêts causales | `DoubleML`, `grf` |

**Modules purement théoriques (13, 14).** Le livrable « implémentation » est
remplacé par des illustrations numériques des bornes / des pathologies. Critère
de passage : dérivations relues et illustrations Monte Carlo cohérentes avec la
théorie (les bornes tiennent, les taux de rejet/couverture correspondent aux
prédictions).

### Extensions économétriques (modules 17–21)

| Fichier | Module | Points clés | Réf. de validation |
|:-------:|--------|-------------|--------------------|
| `17` | **Bootstrap** | plug-in, variance, IC percentile/basique/**BCa**, pairs/résidus | `boot` |
| `18` | **GMM** | moments, GMM efficace, test J, **= 2SLS** | `gmm`, `AER` |
| `19` | **Prédiction conforme** | garantie distribution-libre en échantillon fini | *théorème vérifié* |
| `20` | **Régression quantile** | perte pinball, IRLS, LAD robuste | `quantreg` |
| `21` | **Panel / effets fixes** | within = LSDV (FWL), SE groupées | `plm` |

Ces modules **bouclent des arcs** du projet : le bootstrap complète l'inférence
asymptotique (M1/M2), le GMM unifie IV (M5) et M-estimation (M14), la prédiction
conforme prolonge la validation (M6/M13), et le panel redonne du Frisch-Waugh-Lovell
(M1/M16). Le Module 0 s'enrichit des optimiseurs **Nesterov** et **L-BFGS**.

### Inférence en haute dimension (modules 22–23)

| Fichier | Module | Points clés | Réf. de validation |
|:-------:|--------|-------------|--------------------|
| `22` | **Lasso débiaisé** | projection Zhang-Zhang, score orthogonal nodewise, scaled lasso pour σ, IC valides quand *p > n* | `hdm`, couverture MC |
| `23` | **Analyse de sensibilité (OVB)** | biais de variable omise (Cinelli-Hazlett), R² partiel, **robustness value**, estimation ajustée | `sensemakr` (à 1e-6) |

Le Module 22 **répare** l'inférence post-sélection cassée du Module 14 ; le
Module 23 **quantifie** la fragilité des hypothèses d'identification des Modules
15–16 et 22 (aucune hypothèse remplaçant la randomisation).

### Unification, frontière et fondements (modules 24–28)

| Fichier | Module | Points clés | Réf. de validation |
|:-------:|--------|-------------|--------------------|
| `24` | **Fonctions d'influence** | IC asymptotiquement linéaire ; IC-OLS = sandwich ; = jackknife = bootstrap ; one-step ; EIF de l'ATE + borne semiparamétrique | `sandwich`, M17 |
| `25` | **DiD à adoption échelonnée** | biais du TWFE, poids négatifs (dCDH), Callaway-Sant'Anna, Sun-Abraham | `fixest`, `did` |
| `26` | **Double descente** | seuil d'interpolation, min-norm, seconde descente, ridge monotonise | *phénomène / M0, M4* |
| `27` | **Processus gaussiens & noyaux** | théorème de représentation, GP = kernel ridge, vraisemblance marginale | `DiceKriging` |
| `28` | **Différentiation automatique** | graphe de calcul, mode inverse, backprop (M12) comme cas particulier | `numDeriv` |

Le Module 24 est la **colonne vertébrale** : il montre que M2/M14 (sandwich),
M17 (bootstrap), M16 (Neyman) et M22 (débiaisé) sont une seule idée (la fonction
d'influence). Le Module 25 porte la **frontière économétrique** ; les Modules
26–28 relient sur-paramétrisation, RKHS et la machinerie des frameworks modernes.
La **[synthèse capstone](rapport/capstone.qmd)** tisse les quatre idées qui
traversent les modules.

### Compléments économétriques et statistiques (modules 29–34)

| Fichier | Module | Points clés | Réf. de validation |
|:-------:|--------|-------------|--------------------|
| `29` | **Diagnostics** | hétéroscédasticité (Breusch-Pagan, White), autocorrélation (Durbin-Watson, Breusch-Godfrey), endogénéité (Hausman), suridentification (Sargan), RESET, Jarque-Bera, FGLS | `lmtest`, `AER`, `tseries` |
| `30` | **Variables dép. limitées** | probit, Tobit (censure), Heckman (sélection) | `glm`, `AER`, `sampleSelection` |
| `31` | **Séries temporelles** | ACF/PACF, AR (Yule-Walker), ARMA (CSS), Ljung-Box, Dickey-Fuller | `stats`, `tseries` |
| `32` | **Non paramétrique + RDD** | KDE, Nadaraya-Watson, local linéaire, discontinuité de régression | `stats`, `rdrobust` |
| `33` | **MCMC** | Metropolis-Hastings, Gibbs, Gelman-Rubin, ESS | analytique, `coda` |
| `34` | **Delta + tests multiples** | méthode delta, Bonferroni (FWER), Benjamini-Hochberg (FDR) | `car`, `p.adjust` |

Ces modules **comblent les trous canoniques** d'un socle d'économétrie et de
statistique : détecter (et non seulement corriger) hétéroscédasticité et
endogénéité, traiter les variables limitées, les séries temporelles, le non
paramétrique et l'inférence bayésienne computationnelle. La bibliographie
consolidée est dans [`references.bib`](references.bib) ; un
[glossaire et index des notations](rapport/glossaire.qmd) couvre l'ensemble.

### Applications, études et outils

- **`applications/lalonde.R`** — la boîte à outils causale (OLS, IPW, DML, lasso
  débiaisé, sensibilité) confrontée au **benchmark expérimental NSW** (~1794 \$).
  Reproduit la leçon de LaLonde : sur données observationnelles, la réponse
  dépend de la méthode et du recouvrement.
- **`simulations/study_double_selection.R`** — mini-étude méthodologique
  (Belloni-Chernozhukov-Hansen) : *pourquoi sélectionner les contrôles depuis les
  deux équations ?* La sélection simple s'effondre (couverture 0.38), la double
  reste valide (~0.95).
- **`simulations/benchmark_scratch_vs_ref.R`** — from-scratch vs `lm`/`glm`/
  `ivreg` : précision (1e-15), vitesse, stabilité (QR vs équations normales).
- **`shiny-app/`** — laboratoire interactif (biais-variance, chemins de
  régularisation, orthogonalisation DML) : `shiny::runApp("shiny-app")`.

### Compendium reproductible

Le dépôt est un **compendium de recherche citable** :

- **`renv.lock`** — versions exactes de tous les packages.
- **`_targets.R`** — pipeline à **dépendances suivies** (`make pipeline` /
  `targets::tar_make()`) : seuls les objets impactés sont recalculés ; graphe via
  `targets::tar_visnetwork()`.
- **`Dockerfile`** — environnement bit-à-bit (R épinglé + Quarto + LaTeX, renv
  restauré) : `make docker`.
- **`CITATION.cff`** — métadonnées de citation (prêtes pour un DOI Zenodo).
- **`_pkgdown.yml`** + workflow — **site de documentation** auto-déployé sur
  GitHub Pages (référence des 133 fonctions organisée par thème).
- **`run_tests.R`** — suite complète en processus isolés (cf. ci-dessous).

## Prérequis

- **R ≥ 4.x**
- Packages de validation / outillage : `testthat`, `ggplot2`, `quarto`,
  et les packages de référence par module (`sandwich`, `nlme`, `glmnet`,
  `MASS`, `AER`, `rpart`, `randomForest`, `gbm`, `mclust`, `class`,
  et pour les modules ajoutés `fastshap`/`iml`, `DoubleML`, `grf`).
- **Quarto** pour compiler les dérivations et le rapport.

## Utilisation

```r
# Charger les implémentations
for (f in list.files("R", full.names = TRUE)) source(f)

# Lancer TOUTE la suite de conformité (337 tests) — un processus par fichier
# Rscript run_tests.R

# Reproduire une étude Monte Carlo
source("simulations/mc_ols_gauss_markov.R")
```

> **Exécution des tests.** La suite valide chaque module contre son package de
> référence (`glmnet`, `AER`, `gmm`, `DoubleML`, `grf`, `iml`, `hdm`,
> `sensemakr`…). Chargés *ensemble* dans un même processus R, ces packages lourds
> entrent en conflit (p. ex. `setGeneric` S4 sur `coef`/`mean`, symboles
> exportés, état C accumulé) et corrompent l'environnement — faisant échouer des
> tests **sans rapport** avec le code, alors que chaque fichier passe seul. Le
> script **`run_tests.R`** exécute donc chaque fichier dans un **sous-processus R
> frais** (via `callr`), garantissant l'isolation : `Rscript run_tests.R` →
> *337 pass, 0 fail*. (`testthat::test_dir` en un seul processus reste utilisable
> module par module.)

Compiler une dérivation :

```sh
quarto render derivations/01_ols.qmd
```

## Reproductibilité et qualité

Le projet est structuré en **compendium de recherche reproductible** :

- **Pipeline en une commande** : `make all` (ou `Rscript run_all.R`) charge les
  implémentations, lance les tests et exécute toutes les études Monte Carlo.
  `make book` compile le **livre Quarto** (17 dérivations + rapport + annexe de
  synthèse), `make derivations` / `make rapport` les documents séparément.
- **Environnement verrouillé** : `renv.lock` fige les versions exactes des ~160
  packages (R 4.6.1). `renv::restore()` recrée l'environnement à l'identique.
- **Package documenté** : `NAMESPACE` et 94 pages d'aide `man/*.Rd` générées par
  roxygen2 ; `R CMD check` s'installe et se charge proprement (les avertissements
  résiduels concernent les accents non-ASCII des messages francophones).
- **Rigueur des simulations** : les utilitaires de `R/mc_tools.R` rapportent
  l'**erreur Monte Carlo** de chaque quantité (biais, RMSE, couverture) ;
  `simulations/mc_convergence.R` vérifie la **consistance $\sqrt n$** et les taux.
- **Intégration continue** : `.github/workflows/` fournit `R-CMD-check`
  (Linux + Windows) et la **couverture** (`covr` → Codecov). Couverture ~97 % sur
  les modules noyau.

## Conventions

- R base pour les calculs ; style lisible pour le reste.
- Chaque fonction : docstring roxygen2 (`@param`, `@return`, renvoi à
  l'équation de la dérivation), vérification minimale des entrées.
- `set.seed()` dans toute simulation.
- Commits Git atomiques par livrable (ex. `module 01: dérivation OLS`).

## Références

**Fondamentaux (transverses)**
- Hastie, Tibshirani, Friedman — *The Elements of Statistical Learning* (2ᵉ éd.)
- James, Witten, Hastie, Tibshirani — *An Introduction to Statistical Learning*
- Greene — *Econometric Analysis* (modules 1, 2, 5)
- Boyd & Vandenberghe — *Convex Optimization* (modules 0, 4)
- Bishop — *Pattern Recognition and Machine Learning* (modules 11, 12)

**Théorie de l'apprentissage (module 13)**
- Shalev-Shwartz & Ben-David — *Understanding Machine Learning: From Theory to Algorithms* (référence principale, preuves alignées)
- Mohri, Rostamizadeh & Talwalkar — *Foundations of Machine Learning* (Rademacher)

**Regard économétrique & M-estimation (module 14)**
- van der Vaart — *Asymptotic Statistics* (M-estimation, ch. 5)
- Wooldridge — *Econometric Analysis of Cross Section and Panel Data*
- Berk, Brown, Buja, Zhang, Zhao (2013), « Valid post-selection inference », *Annals of Statistics*

**Interprétabilité (module 15)**
- Lundberg & Lee (2017), « A Unified Approach to Interpreting Model Predictions », NeurIPS
- Molnar — *Interpretable Machine Learning* (en ligne, gratuit)

**Causalité (module 16)**
- Chernozhukov et al. (2018), « Double/debiased machine learning… », *The Econometrics Journal*
- Wager & Athey (2018), « Estimation and Inference of Heterogeneous Treatment Effects using Random Forests », *JASA*
- Athey & Imbens (2016), « Recursive partitioning for heterogeneous causal effects », *PNAS*
