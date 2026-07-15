# Double/Debiased ML pour le modèle partiellement linéaire (éq. 16.4)

Résidualise Y et D par rapport à X (nuisances ML), puis estime
\\\hat\theta = \sum \tilde D\tilde Y / \sum\tilde D^2\\ (score
orthogonal de Neyman). Cross-fitting en K blocs pour débiaiser le
sur-ajustement.

## Usage

``` r
dml_plr(
  y,
  d,
  X,
  K = 5L,
  nuisance = "forest",
  crossfit = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- y:

  réponse.

- d:

  traitement (0/1 ou continu).

- X:

  data.frame (ou matrice) de covariables.

- K:

  nombre de blocs de cross-fitting.

- nuisance:

  "forest" (M9), "boost" (M10) ou "lm".

- crossfit:

  TRUE (cross-fitting) ou FALSE (nuisances sur toutes les données).

- seed:

  graine.

- ...:

  passé aux modèles de nuisance (p.ex. B pour la forêt).

## Value

liste : `theta`, `se`, `ci`, `Ytil`, `Dtil`.
