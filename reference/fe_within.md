# Estimateur à effets fixes (within) pour données de panel (éq. 21.2)

Centre les variables par unité (transformation within, qui élimine
l'effet fixe \\\alpha_i\\) puis applique l'OLS. Renvoie aussi la
variance groupée (clustered) par unité (éq. 21.3). Équivaut au LSDV
(Prop. 21.1).

## Usage

``` r
fe_within(formula, data, id)
```

## Arguments

- formula:

  formule des régresseurs variables (sans les indicatrices).

- data:

  data.frame.

- id:

  nom de la colonne identifiant l'unité.

## Value

liste : `coefficients`, `vcov` (classique within), `vcov_cluster`, `se`,
`se_cluster`, `df.residual`, `residuals`, `N`, `NT`.
