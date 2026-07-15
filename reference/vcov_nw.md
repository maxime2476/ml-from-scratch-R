# Variance HAC de Newey-West (éq. 2.6)

Viande \\\hat\Gamma_0 + \sum\_{\ell=1}^L w\_\ell(\hat\Gamma\_\ell +
\hat\Gamma\_\ell^T)\\ avec poids de Bartlett \\w\_\ell = 1 -
\ell/(L+1)\\. Calé sur
`sandwich::NeweyWest(..., prewhite = FALSE, adjust = FALSE)`.

## Usage

``` r
vcov_nw(fit, lag)
```

## Arguments

- fit:

  objet `ols`.

- lag:

  nombre de retards L.

## Value

matrice de variance p x p.
