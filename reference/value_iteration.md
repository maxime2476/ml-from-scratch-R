# Iteration sur les valeurs (MDP a modele connu) — la reference optimale

Itere l'operateur de Bellman \\Q(s,a)=R(s,a)+\gamma\sum\_{s'}P(s'\|s,a)
\max\_{a'}Q(s',a')\\ jusqu'a convergence. Fournit le \\Q^\\\\ optimal.

## Usage

``` r
value_iteration(P, R, gamma = 0.9, tol = 1e-10)
```

## Arguments

- P:

  transitions, tableau (S x A x S)

- R:

  recompenses (S x A)

- gamma:

  facteur d'actualisation

- tol:

  tolerance.

## Value

liste : `Q`, `V`, `policy`.
