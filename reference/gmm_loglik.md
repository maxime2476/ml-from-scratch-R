# Log-vraisemblance observée d'un mélange gaussien

Log-vraisemblance observée d'un mélange gaussien

## Usage

``` r
gmm_loglik(X, pi, mu, Sigma)
```

## Arguments

- X:

  données n x p.

- pi:

  poids (longueur K).

- mu:

  liste (ou matrice K x p) des moyennes.

- Sigma:

  liste des K covariances.

## Value

la log-vraisemblance observée.
