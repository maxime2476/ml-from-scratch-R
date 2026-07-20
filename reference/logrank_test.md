# Test du log-rank (comparaison de courbes de survie)

Compare les deces **observes** et **attendus** entre groupes a chaque
temps d'evenement ; \\\chi^2\\ sous l'hypothese de survies egales.

## Usage

``` r
logrank_test(time, event, group)
```

## Arguments

- time, event:

  durees et indicateurs

- group:

  facteur a 2 niveaux.

## Value

liste : `statistic`, `df`, `p_value`.
