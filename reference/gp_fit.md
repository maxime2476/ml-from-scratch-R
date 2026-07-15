# Ajustement d'un processus gaussien (régression)

Résout \\(K+\sigma_n^2 I)\alpha=y\\ par Cholesky (Module 0). Renvoie de
quoi prédire moyenne et variance a posteriori, et la **vraisemblance
marginale**.

## Usage

``` r
gp_fit(X, y, lengthscale = 1, sigma_f = 1, sigma_n = 0.1)
```

## Arguments

- X, y:

  données d'apprentissage.

- lengthscale, sigma_f, sigma_n:

  hyperparamètres (échelle, signal, bruit).

## Value

objet `gp` : `alpha`, `L`, `X`, hyperparamètres, `loglik`.
