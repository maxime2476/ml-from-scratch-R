# Complexité de Rademacher empirique d'une classe linéaire à norme bornée

Estime \\\hat{\mathfrak R}\_S = \frac{B}{n}\\\mathbb
E\_\sigma\\\sum_i\sigma_i x_i\\\\ (Déf. 13.10 spécialisée au cas
linéaire) par tirages de signes de Rademacher.

## Usage

``` r
empirical_rademacher_linear(X, B = 1, n_draws = 2000L, seed = NULL)
```

## Arguments

- X:

  matrice n x d des points.

- B:

  rayon de la boule \\\\w\\\_2\le B\\ (défaut 1).

- n_draws:

  nombre de tirages de sigma.

- seed:

  graine.

## Value

l'estimation de la complexité de Rademacher empirique.
