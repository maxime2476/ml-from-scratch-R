# Bagging / forêt aléatoire (bootstrap + agrégation)

Ajuste B arbres (Module 8) sur des rééchantillons bootstrap ; agrège par
moyenne (régression) ou vote majoritaire (classification). Avec `mtry`,
chaque split ne considère qu'un sous-ensemble aléatoire de variables
(forêt aléatoire, décorrélation des arbres). Calcule l'erreur out-of-bag
(éq. 9.2).

## Usage

``` r
bagging_fit(
  formula,
  data,
  method = c("class", "anova"),
  B = 100L,
  mtry = NULL,
  max_depth = 30L,
  min_split = 5L,
  min_leaf = 1L,
  seed = NULL
)
```

## Arguments

- formula:

  formule (prédicteurs numériques).

- data:

  data.frame.

- method:

  "class" ou "anova".

- B:

  nombre d'arbres.

- mtry:

  variables candidates par split ; NULL = bagging (toutes) ; sinon forêt
  aléatoire. Défaut : sqrt(p) (class), p/3 (anova).

- max_depth, min_split, min_leaf:

  hyperparamètres des arbres de base.

- seed:

  graine.

## Value

objet `forest` : `trees`, `oob_idx`, `oob_error`, `oob_pred`, méta.
