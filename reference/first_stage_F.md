# Statistique F de première étape (force des instruments)

Teste la nullité conjointe des coefficients des instruments EXCLUS dans
la régression de première étape du régresseur endogène sur Z (diagnostic
d'instruments faibles ; règle empirique F \< 10, cf. §5 de la
dérivation).

## Usage

``` r
first_stage_F(x_endog, Z, excluded)
```

## Arguments

- x_endog:

  régresseur endogène (vecteur).

- Z:

  matrice des instruments (constante + exogènes + instruments exclus).

- excluded:

  indices des colonnes de Z correspondant aux instruments exclus.

## Value

liste : `F`, `df1`, `df2`, `p_value`.
