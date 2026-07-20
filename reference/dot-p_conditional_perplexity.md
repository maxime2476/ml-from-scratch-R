# Voisinages conditionnels a perplexite fixee (eq. 42.4-42.5)

Calibre un sigma_i PAR POINT pour que chaque ligne de P atteigne la
perplexite cible. L'entropie de P_i croit avec sigma_i (Prop. 42.1),
donc la cible est atteinte en un sigma_i unique, obtenu par dichotomie
sur beta_i = 1/(2 sigma_i^2).

## Usage

``` r
.p_conditional_perplexity(D2, perplexity, tol = 1e-05, max_iter = 50L)
```

## Arguments

- D2:

  matrice n x n des distances au carre.

- perplexity:

  perplexite cible (nombre effectif de voisins).

- tol:

  tolerance sur log2(perplexite).

- max_iter:

  nombre maximal de pas de dichotomie.

## Value

matrice n x n des \\p\_{j\|i}\\ (lignes sommant a 1, diagonale nulle).
