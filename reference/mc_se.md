# Erreur Monte Carlo de la moyenne d'échantillon

\\\mathrm{MCSE}(\bar x) = \hat\sigma/\sqrt R\\. C'est l'incertitude due
au nombre fini R de réplications (à ne pas confondre avec l'écart-type
de l'estimateur étudié).

## Usage

``` r
mc_se(x)
```

## Arguments

- x:

  vecteur des R valeurs simulées.

## Value

l'erreur Monte Carlo de leur moyenne.
