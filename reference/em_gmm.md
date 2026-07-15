# Mélange gaussien par EM (éq. 11.5-11.7)

Étape E : responsabilités (éq. 11.6) ; étape M : mises à jour fermées
(éq. 11.7). Initialisation par k-means. Convergence sur la
log-vraisemblance.

## Usage

``` r
em_gmm(X, K, max_iter = 200L, tol = 1e-08, reg = 1e-06, seed = NULL)
```

## Arguments

- X:

  matrice n x p.

- K:

  nombre de composantes.

- max_iter:

  itérations maximales.

- tol:

  tolérance sur la variation de log-vraisemblance.

- reg:

  régularisation ajoutée à la diagonale des covariances.

- seed:

  graine (initialisation k-means).

## Value

liste : `pi`, `mu` (K x p), `Sigma` (liste), `gamma` (responsabilités),
`loglik`, `iter`, `cluster`.
