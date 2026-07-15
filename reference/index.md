# Package index

## Fondations : algèbre, OLS, inférence

Moindres carrés, factorisations, robustesse, GLM.

- [`solve_ls_chol()`](https://maxime2476.github.io/ml-from-scratch-R/reference/solve_ls_chol.md)
  : Moindres carrés par équations normales et Cholesky
- [`solve_ls_qr()`](https://maxime2476.github.io/ml-from-scratch-R/reference/solve_ls_qr.md)
  : Moindres carrés par QR de Householder
- [`solve_ls_svd()`](https://maxime2476.github.io/ml-from-scratch-R/reference/solve_ls_svd.md)
  : Moindres carrés de norme minimale par SVD
- [`qr_householder()`](https://maxime2476.github.io/ml-from-scratch-R/reference/qr_householder.md)
  : Décomposition QR par réflexions de Householder
- [`ols_confint()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ols_confint.md)
  : Intervalles de confiance des coefficients (éq. 1.8)
- [`ols_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ols_fit.md)
  : Ajustement MCO par QR, avec inférence
- [`ols_ftest()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ols_ftest.md)
  : Test F de restrictions linéaires R beta = r (éq. 1.9)
- [`ols_psi()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ols_psi.md)
  : Fonction d'estimation de l'OLS : psi_i = x_i \* e_i
- [`ols_summary()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ols_summary.md)
  : Tableau récapitulatif d'un ajustement MCO
- [`vcov_hc()`](https://maxime2476.github.io/ml-from-scratch-R/reference/vcov_hc.md)
  : Variances robustes à l'hétéroscédasticité HC0–HC3 (éq. 2.2)
- [`vcov_nw()`](https://maxime2476.github.io/ml-from-scratch-R/reference/vcov_nw.md)
  : Variance HAC de Newey-West (éq. 2.6)
- [`glm_irls()`](https://maxime2476.github.io/ml-from-scratch-R/reference/glm_irls.md)
  : Ajustement d'un GLM par IRLS
- [`back_substitution()`](https://maxime2476.github.io/ml-from-scratch-R/reference/back_substitution.md)
  : Résolution d'un système triangulaire supérieur par remontée
- [`chol_crout()`](https://maxime2476.github.io/ml-from-scratch-R/reference/chol_crout.md)
  : Factorisation de Cholesky (algorithme de Crout)
- [`forward_substitution()`](https://maxime2476.github.io/ml-from-scratch-R/reference/forward_substitution.md)
  : Résolution d'un système triangulaire inférieur par descente
- [`householder_vector()`](https://maxime2476.github.io/ml-from-scratch-R/reference/householder_vector.md)
  : Vecteur de Householder (choix de signe stable)
- [`svd_tools()`](https://maxime2476.github.io/ml-from-scratch-R/reference/svd_tools.md)
  : Outils SVD : rang numérique, conditionnement, pseudo-inverse
- [`coeftest_hc()`](https://maxime2476.github.io/ml-from-scratch-R/reference/coeftest_hc.md)
  : Tests t robustes (coeftest) avec une variance donnée
- [`fwl_beta2()`](https://maxime2476.github.io/ml-from-scratch-R/reference/fwl_beta2.md)
  : Coefficient FWL : régression de M1 y sur M1 X2 (éq. 1.10)
- [`lr_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/lr_test.md)
  : Test du rapport de vraisemblance (éq. 3.12) pour modèles emboîtés
- [`sandwich_vcov()`](https://maxime2476.github.io/ml-from-scratch-R/reference/sandwich_vcov.md)
  : Matrice de variance sandwich générique (éq. 2.1)
- [`score_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/score_test.md)
  : Test du score / Rao (éq. 3.13) pour modèles emboîtés
- [`wald_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/wald_test.md)
  : Test de Wald (éq. 3.11) pour H0 : R beta = r

## Optimisation

Optimiseurs génériques avec preuves de convergence (Module 0).

- [`optim_adam()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_adam.md)
  : Optimiseur Adam (moments adaptatifs)
- [`optim_cd()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_cd.md)
  : Coordinate descent générique
- [`optim_gd()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_gd.md)
  : Descente de gradient à pas constant
- [`optim_lbfgs()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_lbfgs.md)
  : L-BFGS (quasi-Newton à mémoire limitée)
- [`optim_momentum()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_momentum.md)
  : Descente de gradient a momentum (boule pesante de Polyak)
- [`optim_nesterov()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_nesterov.md)
  : Gradient accéléré de Nesterov
- [`optim_newton()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_newton.md)
  : Newton-Raphson (minimisation)
- [`optim_rmsprop()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_rmsprop.md)
  : Optimiseur RMSprop
- [`optim_sgd()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_sgd.md)
  : Gradient stochastique (SGD) par mini-lots

## Régularisation et haute dimension

Ridge, lasso, lasso débiaisé, sensibilité.

- [`adjusted_estimate()`](https://maxime2476.github.io/ml-from-scratch-R/reference/adjusted_estimate.md)
  : Estimation ajustée pour un confondeur de force donnée (OVB)
- [`debiased_lasso()`](https://maxime2476.github.io/ml-from-scratch-R/reference/debiased_lasso.md)
  : Lasso débiaisé / désparsifié (inférence haute dimension valide)
- [`kernel_ridge()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kernel_ridge.md)
  : Régression ridge à noyau (théorème de représentation)
- [`lasso_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/lasso_fit.md)
  : Lasso par coordinate descent (éq. 4.10)
- [`partial_r2()`](https://maxime2476.github.io/ml-from-scratch-R/reference/partial_r2.md)
  : R² partiel d'un effet à partir de sa statistique t
- [`ridge_bias_var()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ridge_bias_var.md)
  : Biais, variance et EQM analytiques du ridge via la SVD (éq. 4.4-4.5)
- [`ridge_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ridge_fit.md)
  : Régression ridge (forme fermée, éq. 4.2)
- [`ridge_posterior()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ridge_posterior.md)
  : Postérieure conjuguée du ridge (Prop. 14.2)
- [`robustness_value()`](https://maxime2476.github.io/ml-from-scratch-R/reference/robustness_value.md)
  : Robustness value (Cinelli-Hazlett 2020)
- [`sensitivity_ols()`](https://maxime2476.github.io/ml-from-scratch-R/reference/sensitivity_ols.md)
  : Analyse de sensibilité complète d'un effet OLS
- [`soft_threshold()`](https://maxime2476.github.io/ml-from-scratch-R/reference/soft_threshold.md)
  : Opérateur de soft-thresholding (éq. 4.9)

## Fonctions d’influence et rééchantillonnage

Le fil unificateur de l’inférence (Modules 17, 24).

- [`influence_mle()`](https://maxime2476.github.io/ml-from-scratch-R/reference/influence_mle.md)
  : Fonction d'influence d'un MLE (GLM canonique)
- [`influence_ols()`](https://maxime2476.github.io/ml-from-scratch-R/reference/influence_ols.md)
  : Fonction d'influence de l'OLS
- [`boot_ci()`](https://maxime2476.github.io/ml-from-scratch-R/reference/boot_ci.md)
  : Intervalles de confiance bootstrap
- [`boot_lm()`](https://maxime2476.github.io/ml-from-scratch-R/reference/boot_lm.md)
  : Bootstrap d'une régression linéaire (pairs ou résidus)
- [`bootstrap()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bootstrap.md)
  : Bootstrap non paramétrique d'une statistique
- [`eif_ate()`](https://maxime2476.github.io/ml-from-scratch-R/reference/eif_ate.md)
  : Fonction d'influence EFFICACE de l'ATE (doublement robuste / AIPW)
- [`jackknife()`](https://maxime2476.github.io/ml-from-scratch-R/reference/jackknife.md)
  : Jackknife (delete-one) — le jackknife infinitésimal EST la fonction
  d'influence
- [`onestep()`](https://maxime2476.github.io/ml-from-scratch-R/reference/onestep.md)
  : Estimateur « en un pas » (one-step / Newton-scoring)

## Causalité

IV, DML, DiD moderne, arbres causaux.

- [`aggregate_att()`](https://maxime2476.github.io/ml-from-scratch-R/reference/aggregate_att.md)
  : Agrégation des ATT(g,t) (Callaway-Sant'Anna)
- [`att_gt()`](https://maxime2476.github.io/ml-from-scratch-R/reference/att_gt.md)
  : ATT groupe-temps (Callaway & Sant'Anna 2021)
- [`causal_tree()`](https://maxime2476.github.io/ml-from-scratch-R/reference/causal_tree.md)
  : Arbre causal minimal (honnête)
- [`dml_plr()`](https://maxime2476.github.io/ml-from-scratch-R/reference/dml_plr.md)
  : Double/Debiased ML pour le modèle partiellement linéaire (éq. 16.4)
- [`first_stage_F()`](https://maxime2476.github.io/ml-from-scratch-R/reference/first_stage_F.md)
  : Statistique F de première étape (force des instruments)
- [`predict_causal_tree()`](https://maxime2476.github.io/ml-from-scratch-R/reference/predict_causal_tree.md)
  : Prédiction du CATE par un arbre causal
- [`sunab()`](https://maxime2476.github.io/ml-from-scratch-R/reference/sunab.md)
  : Event-study de Sun & Abraham (2021) — interaction-weighted
- [`t_learner()`](https://maxime2476.github.io/ml-from-scratch-R/reference/t_learner.md)
  : T-learner : CATE par deux modèles séparés
- [`tsls_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/tsls_fit.md)
  : Estimateur 2SLS (forme fermée, éq. 5.4)
- [`twfe()`](https://maxime2476.github.io/ml-from-scratch-R/reference/twfe.md)
  : Estimateur TWFE (two-way fixed effects) de l'effet du traitement
- [`twfe_weights()`](https://maxime2476.github.io/ml-from-scratch-R/reference/twfe_weights.md)
  : Poids de la régression TWFE (de Chaisemartin & D'Haultfœuille 2020)

## Non linéaire, non supervisé, réseaux

KNN, CART, forêts, boosting, PCA/k-means/EM, MLP, autodiff.

- [`Math(`*`<adnode>`*`)`](https://maxime2476.github.io/ml-from-scratch-R/reference/Math.adnode.md)
  : Fonctions élémentaires enregistrées (exp, log, sqrt, sin, cos, tanh)
- [`Ops(`*`<adnode>`*`)`](https://maxime2476.github.io/ml-from-scratch-R/reference/Ops.adnode.md)
  : Opérateurs élémentaires enregistrés (+, -, \*, /, ^)
- [`Summary(`*`<adnode>`*`)`](https://maxime2476.github.io/ml-from-scratch-R/reference/Summary.adnode.md)
  : Somme enregistrée (réduction scalaire)
- [`ad_cbind1()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ad_cbind1.md)
  : Ajoute une colonne de 1 (terme de biais) — opération enregistrée
- [`ad_grad()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ad_grad.md)
  : Gradient d'une fonction scalaire par mode inverse
- [`ad_reset()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ad_reset.md)
  : Réinitialise la bande d'enregistrement (avant chaque graphe)
- [`adnode()`](https://maxime2476.github.io/ml-from-scratch-R/reference/adnode.md)
  : Crée un nœud de calcul (variable enregistrée sur la bande)
- [`backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/backward.md)
  : Passe arrière : calcule tous les gradients par une seule
  rétropropagation
- [`bagging_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bagging_fit.md)
  : Bagging / forêt aléatoire (bootstrap + agrégation)
- [`batch_norm_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/batch_norm_backward.md)
  : Batch normalization (passe arriere)
- [`best_split()`](https://maxime2476.github.io/ml-from-scratch-R/reference/best_split.md)
  : Meilleur split d'un nœud (éq. 8.3)
- [`boost_loss_path()`](https://maxime2476.github.io/ml-from-scratch-R/reference/boost_loss_path.md)
  : Trajectoire de la perte d'entraînement/test selon le nombre d'arbres
- [`cart_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/cart_fit.md)
  : Ajustement d'un arbre CART
- [`conv2d_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/conv2d_backward.md)
  : Retropropagation de la convolution 2D
- [`cost_complexity_prune()`](https://maxime2476.github.io/ml-from-scratch-R/reference/cost_complexity_prune.md)
  : Élagage coût-complexité (éq. 8.4-8.5)
- [`em_gmm()`](https://maxime2476.github.io/ml-from-scratch-R/reference/em_gmm.md)
  : Mélange gaussien par EM (éq. 11.5-11.7)
- [`gradient_boost()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gradient_boost.md)
  : Gradient boosting (descente de gradient fonctionnelle, éq.
  10.1-10.2)
- [`impurity_entropy()`](https://maxime2476.github.io/ml-from-scratch-R/reference/impurity_entropy.md)
  : Entropie d'un vecteur d'étiquettes (éq. 8.1)
- [`impurity_gini()`](https://maxime2476.github.io/ml-from-scratch-R/reference/impurity_gini.md)
  : Impureté de Gini d'un vecteur d'étiquettes (éq. 8.1)
- [`impurity_variance()`](https://maxime2476.github.io/ml-from-scratch-R/reference/impurity_variance.md)
  : Impureté de variance (régression, éq. 8.2)
- [`kernel_pca()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kernel_pca.md)
  : Analyse en composantes principales a noyau (kernel PCA)
- [`kmeans_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kmeans_fit.md)
  : k-means par l'algorithme de Lloyd (éq. 11.4)
- [`knn_classify()`](https://maxime2476.github.io/ml-from-scratch-R/reference/knn_classify.md)
  : Classification KNN par vote majoritaire
- [`knn_regression()`](https://maxime2476.github.io/ml-from-scratch-R/reference/knn_regression.md)
  : Régression KNN (éq. 7.1)
- [`max_pool2d_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/max_pool2d_backward.md)
  : Retropropagation du max-pooling
- [`mlp_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mlp_backward.md)
  : Rétropropagation : gradients analytiques (éq. 12.2-12.5)
- [`mlp_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mlp_fit.md)
  : Entraînement d'un MLP par SGD (éq. 12.6)
- [`mlp_forward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mlp_forward.md)
  : Passe avant du MLP (éq. 12.1)
- [`mlp_loss()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mlp_loss.md)
  : Perte moyenne du MLP
- [`mlp_numgrad()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mlp_numgrad.md)
  : Gradient numérique par différences finies centrées (éq. 12.7)
- [`mm()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mm.md)
  : Produit matriciel enregistré
- [`multi_head_attention()`](https://maxime2476.github.io/ml-from-scratch-R/reference/multi_head_attention.md)
  : Attention multi-tetes
- [`n_leaves()`](https://maxime2476.github.io/ml-from-scratch-R/reference/n_leaves.md)
  : Nombre de feuilles d'un arbre CART
- [`pca_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/pca_fit.md)
  : Analyse en composantes principales par SVD (éq. 11.2-11.3)
- [`predict_boost()`](https://maxime2476.github.io/ml-from-scratch-R/reference/predict_boost.md)
  : Prédiction d'un modèle de boosting
- [`predict_cart()`](https://maxime2476.github.io/ml-from-scratch-R/reference/predict_cart.md)
  : Prédiction d'un arbre CART
- [`predict_forest()`](https://maxime2476.github.io/ml-from-scratch-R/reference/predict_forest.md)
  : Prédiction d'une forêt / d'un ensemble baggé
- [`predict_mlp()`](https://maxime2476.github.io/ml-from-scratch-R/reference/predict_mlp.md)
  : Prédiction d'un MLP
- [`random_forest_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/random_forest_fit.md)
  : Forêt aléatoire (alias de bagging_fit avec mtry actif)
- [`rnn_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/rnn_backward.md)
  : Retropropagation dans le temps (BPTT) du RNN simple

## Noyaux et processus gaussiens

- [`double_descent_curve()`](https://maxime2476.github.io/ml-from-scratch-R/reference/double_descent_curve.md)
  : Courbe de risque en fonction de la complexité (double descente)
- [`fit_rff()`](https://maxime2476.github.io/ml-from-scratch-R/reference/fit_rff.md)
  : Ajuste un modèle à caractéristiques aléatoires (interpolant ou
  régularisé)
- [`gp_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gp_fit.md)
  : Ajustement d'un processus gaussien (régression)
- [`gp_optimize()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gp_optimize.md)
  : Sélection des hyperparamètres par maximum de vraisemblance marginale
- [`gp_predict()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gp_predict.md)
  : Prédiction d'un processus gaussien
- [`kernel_ridge()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kernel_ridge.md)
  : Régression ridge à noyau (théorème de représentation)
- [`random_features()`](https://maxime2476.github.io/ml-from-scratch-R/reference/random_features.md)
  : Caractéristiques de Fourier aléatoires (approximation du noyau
  gaussien)
- [`rbf_kernel()`](https://maxime2476.github.io/ml-from-scratch-R/reference/rbf_kernel.md)
  : Noyau gaussien (RBF)

## Interprétabilité et validation

- [`bw_loocv()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bw_loocv.md)
  : Selection de fenetre par validation croisee leave-one-out
  (Nadaraya-Watson)
- [`conformal_quantile()`](https://maxime2476.github.io/ml-from-scratch-R/reference/conformal_quantile.md)
  : Quantile conforme (éq. 19.1)
- [`conformal_split()`](https://maxime2476.github.io/ml-from-scratch-R/reference/conformal_split.md)
  : Prédiction conforme par découpage (split conformal)
- [`gcv_linear()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gcv_linear.md)
  : Validation croisée généralisée (GCV)
- [`ice()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ice.md)
  : Espérances conditionnelles individuelles (ICE, éq. 15.3)
- [`info_criteria()`](https://maxime2476.github.io/ml-from-scratch-R/reference/info_criteria.md)
  : AIC et BIC d'un modèle ajusté (éq. 6.5, 6.7)
- [`kfold_cv()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kfold_cv.md)
  : Validation croisée K-fold générique (éq. 6.2)
- [`loocv_linear()`](https://maxime2476.github.io/ml-from-scratch-R/reference/loocv_linear.md)
  : LOOCV fermé pour la régression linéaire (éq. 6.3)
- [`mallows_cp()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mallows_cp.md)
  : Cp de Mallows (éq. 6.4)
- [`pdp()`](https://maxime2476.github.io/ml-from-scratch-R/reference/pdp.md)
  : Dépendance partielle (PDP, éq. 15.2)
- [`permutation_importance()`](https://maxime2476.github.io/ml-from-scratch-R/reference/permutation_importance.md)
  : Importance par permutation (éq. 15.6)
- [`pinball_loss()`](https://maxime2476.github.io/ml-from-scratch-R/reference/pinball_loss.md)
  : Perte pinball (fonction « check », éq. 20.1)
- [`qreg_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/qreg_fit.md)
  : Régression quantile par IRLS (éq. 20.2)
- [`shapley_exact()`](https://maxime2476.github.io/ml-from-scratch-R/reference/shapley_exact.md)
  : Valeurs de Shapley exactes par énumération (éq. 15.4, p \<= 10)
- [`shapley_permutation()`](https://maxime2476.github.io/ml-from-scratch-R/reference/shapley_permutation.md)
  : Valeurs de Shapley approchées par échantillonnage de permutations

## Économétrie : GMM, panel, quantile

- [`em_gmm()`](https://maxime2476.github.io/ml-from-scratch-R/reference/em_gmm.md)
  : Mélange gaussien par EM (éq. 11.5-11.7)
- [`fe_within()`](https://maxime2476.github.io/ml-from-scratch-R/reference/fe_within.md)
  : Estimateur à effets fixes (within) pour données de panel (éq. 21.2)
- [`gls_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gls_fit.md)
  : Moindres carrés généralisés (GLS, éq. 2.4) avec Omega connue
- [`gmm_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gmm_fit.md)
  : GMM générique (moments non linéaires) par optimisation
- [`gmm_linear()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gmm_linear.md)
  : GMM linéaire (moments d'instruments) : 2SLS et GMM efficace à deux
  étapes
- [`gmm_loglik()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gmm_loglik.md)
  : Log-vraisemblance observée d'un mélange gaussien
- [`wls_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/wls_fit.md)
  : Moindres carrés pondérés (WLS, éq. 2.5)

## Outils Monte Carlo et théorie de l’apprentissage

- [`bias_variance_mc()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bias_variance_mc.md)
  : Estimation Monte Carlo de la décomposition biais-variance (éq. 6.1)
- [`convergence_study()`](https://maxime2476.github.io/ml-from-scratch-R/reference/convergence_study.md)
  : Étude de convergence : biais/RMSE et diagnostics de taux selon n
- [`coverage_mc()`](https://maxime2476.github.io/ml-from-scratch-R/reference/coverage_mc.md)
  : Couverture empirique d'IC avec son erreur Monte Carlo (binomiale)
- [`distance_concentration()`](https://maxime2476.github.io/ml-from-scratch-R/reference/distance_concentration.md)
  : Concentration des distances (éq. 7.4-7.5)
- [`edge_length()`](https://maxime2476.github.io/ml-from-scratch-R/reference/edge_length.md)
  : Longueur d'arête d'un voisinage capturant une fraction r (éq. 7.3)
- [`empirical_rademacher_linear()`](https://maxime2476.github.io/ml-from-scratch-R/reference/empirical_rademacher_linear.md)
  : Complexité de Rademacher empirique d'une classe linéaire à norme
  bornée
- [`gaussian_loglik()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gaussian_loglik.md)
  : Log-vraisemblance gaussienne (régression linéaire) au MLE
- [`hoeffding_bound()`](https://maxime2476.github.io/ml-from-scratch-R/reference/hoeffding_bound.md)
  : Borne de Hoeffding bilatérale (éq. 13.2)
- [`is_separable()`](https://maxime2476.github.io/ml-from-scratch-R/reference/is_separable.md)
  : Séparabilité linéaire d'un étiquetage (théorème de Gordan)
- [`m_estimation_vcov()`](https://maxime2476.github.io/ml-from-scratch-R/reference/m_estimation_vcov.md)
  : Variance sandwich d'un M-estimateur (éq. 14.2)
- [`mc_se()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mc_se.md)
  : Erreur Monte Carlo de la moyenne d'échantillon
- [`mc_summary()`](https://maxime2476.github.io/ml-from-scratch-R/reference/mc_summary.md)
  : Résumé d'une étude de simulation : biais, RMSE, variance, avec
  erreurs MC
- [`ols_psi()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ols_psi.md)
  : Fonction d'estimation de l'OLS : psi_i = x_i \* e_i
- [`rademacher_linear_bound()`](https://maxime2476.github.io/ml-from-scratch-R/reference/rademacher_linear_bound.md)
  : Borne théorique de Rademacher pour une classe linéaire (éq. 13.5)
- [`reject_mc()`](https://maxime2476.github.io/ml-from-scratch-R/reference/reject_mc.md)
  : Taux de rejet (taille/puissance) avec erreur Monte Carlo
- [`ridge_posterior()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ridge_posterior.md)
  : Postérieure conjuguée du ridge (Prop. 14.2)
- [`rmse_rate()`](https://maxime2476.github.io/ml-from-scratch-R/reference/rmse_rate.md)
  : Pente log-log du RMSE en fonction de n (taux de convergence
  empirique)
- [`shatters_hyperplane()`](https://maxime2476.github.io/ml-from-scratch-R/reference/shatters_hyperplane.md)
  : La classe des hyperplans pulvérise-t-elle un ensemble de points ?
  (Déf. 13.6)

## Diagnostics, tests et outils d’inférence

Tests de spécification (Module 29), méthode delta et tests multiples
(Module 34).

- [`bg_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bg_test.md)
  : Test de Breusch-Godfrey (autocorrelation d'ordre p, LM)
- [`bp_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bp_test.md)
  : Test de Breusch-Pagan (heteroscedasticite)
- [`delta_method()`](https://maxime2476.github.io/ml-from-scratch-R/reference/delta_method.md)
  : Methode delta : variance asymptotique d'une fonction d'un estimateur
- [`dw_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/dw_test.md)
  : Statistique de Durbin-Watson (autocorrelation d'ordre 1)
- [`dwh_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/dwh_test.md)
  : Test d'endogeneite de Durbin-Wu-Hausman (regression augmentee)
- [`fgls()`](https://maxime2476.github.io/ml-from-scratch-R/reference/fgls.md)
  : Moindres carres generalises FAISABLES (FGLS) pour heteroscedasticite
- [`jarque_bera()`](https://maxime2476.github.io/ml-from-scratch-R/reference/jarque_bera.md)
  : Test de normalite de Jarque-Bera
- [`p_adjust_bh()`](https://maxime2476.github.io/ml-from-scratch-R/reference/p_adjust_bh.md)
  : Correction de Benjamini-Hochberg (controle du FDR)
- [`p_adjust_bonferroni()`](https://maxime2476.github.io/ml-from-scratch-R/reference/p_adjust_bonferroni.md)
  : Correction de Bonferroni (controle du FWER)
- [`reset_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/reset_test.md)
  : Test RESET de Ramsey (forme fonctionnelle)
- [`sargan_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/sargan_test.md)
  : Test de suridentification de Sargan
- [`white_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/white_test.md)
  : Test de White (heteroscedasticite generale)

## Variables dépendantes limitées

- [`heckman()`](https://maxime2476.github.io/ml-from-scratch-R/reference/heckman.md)
  : Modele de selection de Heckman (estimation en deux etapes)
- [`probit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/probit.md)
  : Probit (reponse binaire) par IRLS (lien probit)
- [`tobit_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/tobit_fit.md)
  : Regression Tobit (reponse censuree) par maximum de vraisemblance

## Séries temporelles

- [`acf_ts()`](https://maxime2476.github.io/ml-from-scratch-R/reference/acf_ts.md)
  : Fonction d'autocorrelation (ACF)
- [`adf_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/adf_test.md)
  : Test de racine unitaire de Dickey-Fuller augmente (ADF)
- [`ar_yw()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ar_yw.md)
  : Estimation d'un AR(p) par les equations de Yule-Walker
- [`arma_css()`](https://maxime2476.github.io/ml-from-scratch-R/reference/arma_css.md)
  : Estimation d'un ARMA(p,q) par moindres carres conditionnels (CSS)
- [`ljung_box()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ljung_box.md)
  : Test portmanteau de Ljung-Box (autocorrelation residuelle)
- [`pacf_ts()`](https://maxime2476.github.io/ml-from-scratch-R/reference/pacf_ts.md)
  : Fonction d'autocorrelation partielle (PACF, recursion de
  Durbin-Levinson)

## Non paramétrique, RDD et MCMC

- [`bw_loocv()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bw_loocv.md)
  : Selection de fenetre par validation croisee leave-one-out
  (Nadaraya-Watson)
- [`ess()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ess.md)
  : Taille d'echantillon effective (ESS)
- [`gelman_rubin()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gelman_rubin.md)
  : Diagnostic de convergence de Gelman-Rubin (statistique R-hat)
- [`gibbs_linreg()`](https://maxime2476.github.io/ml-from-scratch-R/reference/gibbs_linreg.md)
  : Gibbs pour la regression lineaire bayesienne (prior non informatif)
- [`kde()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kde.md)
  : Estimation de densite par noyau (KDE)
- [`local_linear()`](https://maxime2476.github.io/ml-from-scratch-R/reference/local_linear.md)
  : Regression locale lineaire
- [`metropolis_hastings()`](https://maxime2476.github.io/ml-from-scratch-R/reference/metropolis_hastings.md)
  : Echantillonneur de Metropolis-Hastings (marche aleatoire)
- [`nadaraya_watson()`](https://maxime2476.github.io/ml-from-scratch-R/reference/nadaraya_watson.md)
  : Regression de Nadaraya-Watson (local constant)
- [`rdd()`](https://maxime2476.github.io/ml-from-scratch-R/reference/rdd.md)
  : Discontinuite de regression (RDD, effet local au seuil)

## Réseaux de neurones profonds

Optimiseurs, régularisation, CNN, RNN/LSTM, attention (Modules 35-38).

- [`attention()`](https://maxime2476.github.io/ml-from-scratch-R/reference/attention.md)
  : Attention par produit scalaire mis a l'echelle
- [`batch_norm()`](https://maxime2476.github.io/ml-from-scratch-R/reference/batch_norm.md)
  : Batch normalization (passe avant)
- [`batch_norm_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/batch_norm_backward.md)
  : Batch normalization (passe arriere)
- [`conv2d()`](https://maxime2476.github.io/ml-from-scratch-R/reference/conv2d.md)
  : Convolution 2D (correlation croisee, mode "valid")
- [`conv2d_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/conv2d_backward.md)
  : Retropropagation de la convolution 2D
- [`dropout()`](https://maxime2476.github.io/ml-from-scratch-R/reference/dropout.md)
  : Dropout (inverted dropout)
- [`layer_norm()`](https://maxime2476.github.io/ml-from-scratch-R/reference/layer_norm.md)
  : Normalisation de couche (layer norm)
- [`lstm_forward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/lstm_forward.md)
  : Passe avant d'une cellule LSTM (sequence)
- [`max_pool2d()`](https://maxime2476.github.io/ml-from-scratch-R/reference/max_pool2d.md)
  : Max-pooling 2D (sous-echantillonnage par le maximum)
- [`max_pool2d_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/max_pool2d_backward.md)
  : Retropropagation du max-pooling
- [`multi_head_attention()`](https://maxime2476.github.io/ml-from-scratch-R/reference/multi_head_attention.md)
  : Attention multi-tetes
- [`optim_adam()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_adam.md)
  : Optimiseur Adam (moments adaptatifs)
- [`optim_momentum()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_momentum.md)
  : Descente de gradient a momentum (boule pesante de Polyak)
- [`optim_rmsprop()`](https://maxime2476.github.io/ml-from-scratch-R/reference/optim_rmsprop.md)
  : Optimiseur RMSprop
- [`positional_encoding()`](https://maxime2476.github.io/ml-from-scratch-R/reference/positional_encoding.md)
  : Encodage positionnel sinusoidal
- [`rnn_backward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/rnn_backward.md)
  : Retropropagation dans le temps (BPTT) du RNN simple
- [`rnn_forward()`](https://maxime2476.github.io/ml-from-scratch-R/reference/rnn_forward.md)
  : Passe avant d'un RNN simple (sequence -\> sequence)
- [`softmax_rows()`](https://maxime2476.github.io/ml-from-scratch-R/reference/softmax_rows.md)
  : Softmax numeriquement stable (par ligne)

## Machine learning : SVM, génératifs, clustering, réduction de dimension, RL

Modules 39-43.

- [`agglomerative()`](https://maxime2476.github.io/ml-from-scratch-R/reference/agglomerative.md)
  : Clustering hierarchique agglomeratif (assignation a k groupes)
- [`bandit_thompson()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bandit_thompson.md)
  : Bandit par echantillonnage de Thompson (Bernoulli)
- [`bandit_ucb()`](https://maxime2476.github.io/ml-from-scratch-R/reference/bandit_ucb.md)
  : Bandit UCB1 (borne de confiance superieure)
- [`dbscan_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/dbscan_fit.md)
  : DBSCAN (clustering base sur la densite)
- [`ica_fastica()`](https://maxime2476.github.io/ml-from-scratch-R/reference/ica_fastica.md)
  : Analyse en composantes independantes (FastICA)
- [`kernel_pca()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kernel_pca.md)
  : Analyse en composantes principales a noyau (kernel PCA)
- [`lda_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/lda_fit.md)
  : Analyse discriminante lineaire (LDA)
- [`lda_predict()`](https://maxime2476.github.io/ml-from-scratch-R/reference/lda_predict.md)
  : Prediction LDA
- [`naive_bayes_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/naive_bayes_fit.md)
  : Classifieur naif de Bayes (gaussien)
- [`naive_bayes_predict()`](https://maxime2476.github.io/ml-from-scratch-R/reference/naive_bayes_predict.md)
  : Prediction Naive Bayes
- [`nmf()`](https://maxime2476.github.io/ml-from-scratch-R/reference/nmf.md)
  : Factorisation en matrices non negatives (NMF)
- [`q_learning()`](https://maxime2476.github.io/ml-from-scratch-R/reference/q_learning.md)
  : Q-learning tabulaire (MDP a modele INCONNU)
- [`qda_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/qda_fit.md)
  : Analyse discriminante quadratique (QDA)
- [`qda_predict()`](https://maxime2476.github.io/ml-from-scratch-R/reference/qda_predict.md)
  : Prediction QDA
- [`spectral_clustering()`](https://maxime2476.github.io/ml-from-scratch-R/reference/spectral_clustering.md)
  : Clustering spectral
- [`svm_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/svm_fit.md)
  : SVM a marge souple par le probleme dual (programmation quadratique)
- [`svm_predict()`](https://maxime2476.github.io/ml-from-scratch-R/reference/svm_predict.md)
  : Prediction d'un SVM
- [`svm_rbf()`](https://maxime2476.github.io/ml-from-scratch-R/reference/svm_rbf.md)
  : Noyau gaussien (RBF) pour SVM
- [`tsne()`](https://maxime2476.github.io/ml-from-scratch-R/reference/tsne.md)
  : t-SNE (visualisation preservant le voisinage) — version compacte
- [`value_iteration()`](https://maxime2476.github.io/ml-from-scratch-R/reference/value_iteration.md)
  : Iteration sur les valeurs (MDP a modele connu) — la reference
  optimale

## Économétrie avancée : VAR, GARCH, survie, panel

Modules 44-47.

- [`arch_lm_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/arch_lm_test.md)
  : Test ARCH-LM (heteroscedasticite conditionnelle)
- [`cox_ph()`](https://maxime2476.github.io/ml-from-scratch-R/reference/cox_ph.md)
  : Modele de Cox a risques proportionnels (vraisemblance partielle)
- [`dynamic_panel_fe()`](https://maxime2476.github.io/ml-from-scratch-R/reference/dynamic_panel_fe.md)
  : Estimateur within (effets fixes) d'un panel dynamique — pour
  illustrer le biais
- [`dynamic_panel_iv()`](https://maxime2476.github.io/ml-from-scratch-R/reference/dynamic_panel_iv.md)
  : Panel dynamique par variables instrumentales (Anderson-Hsiao /
  Arellano-Bond)
- [`garch_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/garch_fit.md)
  : Estimation d'un GARCH(1,1) par (quasi-)maximum de vraisemblance
- [`granger_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/granger_test.md)
  : Test de causalite de Granger
- [`kaplan_meier()`](https://maxime2476.github.io/ml-from-scratch-R/reference/kaplan_meier.md)
  : Estimateur de survie de Kaplan-Meier
- [`logrank_test()`](https://maxime2476.github.io/ml-from-scratch-R/reference/logrank_test.md)
  : Test du log-rank (comparaison de courbes de survie)
- [`synthetic_control()`](https://maxime2476.github.io/ml-from-scratch-R/reference/synthetic_control.md)
  : Controle synthetique (Abadie)
- [`var_fit()`](https://maxime2476.github.io/ml-from-scratch-R/reference/var_fit.md)
  : Estimation d'un VAR(p) (equation par equation, OLS)
- [`var_irf()`](https://maxime2476.github.io/ml-from-scratch-R/reference/var_irf.md)
  : Fonctions de reponse impulsionnelle (IRF) orthogonalisees

## Package
