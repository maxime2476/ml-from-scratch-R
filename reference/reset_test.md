# Test RESET de Ramsey (forme fonctionnelle)

Ajoute des puissances des valeurs ajustees \\\hat y^2,\hat y^3,\dots\\
au modele et teste (F) leur significativite jointe. Rejet =
non-linearite negligee.

## Usage

``` r
reset_test(formula, data, powers = 2:3)
```

## Arguments

- formula, data:

  cf. `bp_test`

- powers:

  puissances de \\\hat y\\

## Value

liste : `statistic` (F), `df1`, `df2`, `p_value`.
