# Taux de rejet (taille/puissance) avec erreur Monte Carlo

Taux de rejet (taille/puissance) avec erreur Monte Carlo

## Usage

``` r
reject_mc(rejected, nominal = 0.05)
```

## Arguments

- rejected:

  vecteur logique des rejets.

- nominal:

  niveau nominal du test (défaut 0.05, pour la taille).

## Value

liste : `rate`, `se`, `ci`, `R`, `nominal_ok`.
