# Gradient numérique par différences finies centrées (éq. 12.7)

Sert à VÉRIFIER la rétropropagation. Perturbe chaque paramètre de +/-
eps et approche le gradient par
\\(L(\theta+\varepsilon)-L(\theta-\varepsilon))/(2\varepsilon)\\.

## Usage

``` r
mlp_numgrad(params, X, Y, activation, loss, eps = 1e-06)
```

## Arguments

- params:

  liste des paramètres.

- X, Y:

  données.

- activation:

  activation cachée.

- loss:

  "mse" ou "logloss".

- eps:

  pas de différence finie.

## Value

liste des gradients numériques (même structure que `params`).
