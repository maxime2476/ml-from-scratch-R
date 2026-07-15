# Lasso débiaisé / désparsifié (inférence haute dimension valide)

Corrige le biais de rétrécissement du lasso par une **projection de
faible dimension** (Zhang-Zhang 2014) : pour chaque coordonnée j, on
partiale les autres covariables hors de \\x_j\\ par un lasso « nodewise
», et on utilise le résidu \\\hat\tau_j\\ comme score orthogonal (même
idée que Neyman, M16) : \\\hat\beta^d_j=\hat\beta_j +
\hat\tau_j^\top(y-X\hat\beta)/(\hat\tau_j^\top x_j)\\. Fournit des
intervalles de confiance **valides** même si p \> n — là où le t-test
post-lasso naïf échoue (Module 14).

## Usage

``` r
debiased_lasso(
  X,
  y,
  lambda = NULL,
  lambda_node = NULL,
  targets = NULL,
  sigma = NULL,
  level = 0.95
)
```

## Arguments

- X:

  matrice n x p.

- y:

  réponse.

- lambda:

  pénalité du lasso principal (défaut : théorique).

- lambda_node:

  pénalité des lasso nodewise (défaut : théorique).

- targets:

  indices des coordonnées à débiaiser (défaut : toutes).

- sigma:

  écart-type du bruit (défaut : estimé sur les résidus du lasso).

- level:

  niveau de confiance des intervalles (défaut 0.95).

## Value

liste : `estimate`, `se`, `lower`, `upper`, `beta_lasso`, `sigma`.
