# Gradient boosting (descente de gradient fonctionnelle, éq. 10.1-10.2)

Ajuste additivement M arbres aux pseudo-résidus. `loss = "l2"` (résidu,
éq. 10.3) ou `"logloss"` (y - p, éq. 10.5). Avec `newton = TRUE`, les
poids de feuille sont ceux de Newton (éq. 10.7) — pour la log-loss, la
mise à jour de Friedman (comme `gbm`).

## Usage

``` r
gradient_boost(
  formula,
  data,
  loss = c("l2", "logloss"),
  M = 100L,
  nu = 0.1,
  max_depth = 3L,
  min_split = 10L,
  min_leaf = 5L,
  lambda = 0,
  newton = TRUE
)
```

## Arguments

- formula:

  formule (prédicteurs numériques).

- data:

  data.frame.

- loss:

  "l2" (régression) ou "logloss" (classification binaire `y in {0,1}`).

- M:

  nombre d'arbres.

- nu:

  taux d'apprentissage (shrinkage) dans (0,1\].

- max_depth, min_split, min_leaf:

  hyperparamètres des arbres de base.

- lambda:

  régularisation L2 des poids de feuille (éq. 10.7).

- newton:

  raffiner les feuilles par le pas de Newton (défaut TRUE).

## Value

objet `boost`.
