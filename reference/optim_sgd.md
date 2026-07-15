# Gradient stochastique (SGD) par mini-lots

Minimise une perte-somme \\f(x)=\frac1n\sum_i f_i(x)\\ par l'itération
\\x\_{k+1}=x_k - t_k g_k\\ avec \\g_k\\ gradient sur un mini-lot (éq.
0.20). Pas constant (`step`) ou décroissant (`step_fun`, conditions de
Robbins-Monro 0.21).

## Usage

``` r
optim_sgd(
  grad_i,
  x0,
  n,
  batch = 1L,
  step = NULL,
  step_fun = NULL,
  epochs = 50L,
  seed = NULL
)
```

## Arguments

- grad_i:

  fonction : `grad_i(x, idx)` renvoie le gradient moyen sur les
  observations d'indices `idx`.

- x0:

  point initial.

- n:

  nombre total d'observations.

- batch:

  taille de mini-lot.

- step:

  pas constant (ignoré si `step_fun` fourni).

- step_fun:

  (optionnel) fonction du compteur d'updates t -\> pas.

- epochs:

  nombre de passages sur l'échantillon.

- seed:

  (optionnel) graine pour la permutation des indices.

## Value

liste : `par`, `updates`.
