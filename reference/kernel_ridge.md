# Régression ridge à noyau (théorème de représentation)

\\\hat f(x)=k(x)^\top(K+\lambda I)^{-1}y\\. **Identique** à la moyenne a
posteriori du GP avec \\\lambda=\sigma_n^2\\ (pont
bayésien/fréquentiste).

## Usage

``` r
kernel_ridge(X, y, lengthscale = 1, sigma_f = 1, lambda = 0.01)
```

## Arguments

- X, y:

  données.

- lengthscale, sigma_f:

  paramètres du noyau.

- lambda:

  pénalité.

## Value

fonction `newX -> prédictions`.
