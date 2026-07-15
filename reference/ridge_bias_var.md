# Biais, variance et EQM analytiques du ridge via la SVD (éq. 4.4-4.5)

Décompose l'EQM du ridge en biais² + variance dans la base des
composantes principales (rotation V de la SVD de X). Sert à confronter
la théorie au Monte Carlo.

## Usage

``` r
ridge_bias_var(X, beta_true, sigma2, lambda)
```

## Arguments

- X:

  design (standardisé si l'on veut la cohérence avec `ridge_fit`).

- beta_true:

  vecteur des coefficients vrais (échelle des colonnes de X).

- sigma2:

  variance des erreurs.

- lambda:

  pénalité.

## Value

liste : `bias2`, `variance`, `mse` (totaux), et vecteurs par composante.
