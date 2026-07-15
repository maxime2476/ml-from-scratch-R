# Étude de convergence : biais/RMSE et diagnostics de taux selon n

Pour chaque taille n, exécute R réplications de `sim_fn(n)` (qui renvoie
un \\\hat\theta\\) et calcule biais, RMSE (avec erreurs MC), ainsi que
les diagnostics de **taux** \\\sqrt n\\\text{biais}\\ (doit rester borné
si le biais est \\O(1/\sqrt n)\\ ou mieux) et \\\sqrt n\\\hat{sd}\\
(doit se stabiliser vers l'écart-type asymptotique si l'estimateur est
\\\sqrt n\\- consistant).

## Usage

``` r
convergence_study(sim_fn, ns, R, truth, seed = NULL)
```

## Arguments

- sim_fn:

  fonction `n -> theta_hat` (une réplication).

- ns:

  vecteur des tailles d'échantillon.

- R:

  réplications par taille.

- truth:

  valeur vraie.

- seed:

  graine.

## Value

data.frame : `n`, `bias`, `bias_se`, `rmse`, `rmse_se`, `sqrtn_bias`,
`sqrtn_sd`.
