# Couverture empirique d'IC avec son erreur Monte Carlo (binomiale)

Une couverture estimée sur R réplications est une proportion : son
erreur MC est \\\sqrt{\hat p(1-\hat p)/R}\\. On peut ainsi juger si
l'écart à la valeur nominale est significatif.

## Usage

``` r
coverage_mc(covered, nominal = 0.95)
```

## Arguments

- covered:

  vecteur logique (l'IC contenait la vraie valeur).

- nominal:

  niveau nominal (défaut 0.95).

## Value

liste : `coverage`, `se`, `ci` (IC de la couverture), `R`, `nominal`,
`nominal_ok` (TRUE si l'IC contient le niveau nominal).
