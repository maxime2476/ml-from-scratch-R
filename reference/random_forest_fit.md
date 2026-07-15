# Forêt aléatoire (alias de bagging_fit avec mtry actif)

Forêt aléatoire (alias de bagging_fit avec mtry actif)

## Usage

``` r
random_forest_fit(
  formula,
  data,
  method = c("class", "anova"),
  B = 100L,
  mtry = NULL,
  ...
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

- ...:

  arguments supplémentaires transmis à `bagging_fit`.

## Value

objet `forest`.
