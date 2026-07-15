# Estimation d'un GARCH(1,1) par (quasi-)maximum de vraisemblance

\\x_t=\sigma_t z_t\\, \\z_t\sim(0,1)\\, avec la variance conditionnelle
\\\sigma_t^2=\omega+\alpha x\_{t-1}^2+\beta\sigma\_{t-1}^2\\. On
maximise la log-vraisemblance gaussienne (paramétrage assurant
\\\omega\>0\\, \\\alpha, \beta\ge0\\, \\\alpha+\beta\<1\\ :
stationnarite).

## Usage

``` r
garch_fit(x, maxit = 500L)
```

## Arguments

- x:

  serie (rendements centres) ; @param maxit iterations de l'optimiseur.

## Value

liste : `omega`, `alpha`, `beta`, `sigma` (volatilite conditionnelle),
`loglik`.
