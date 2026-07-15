# Régression quantile par IRLS (éq. 20.2)

Minimise la perte pinball par moindres carrés pondérés itérés
(majoration de Hunter-Lange). Réutilise la QR pondérée du Module 0.
Reproche [`quantreg::rq`](https://rdrr.io/pkg/quantreg/man/rq.html).

## Usage

``` r
qreg_fit(formula, data, tau = 0.5, maxit = 200L, tol = 1e-08, eps = 1e-06)
```

## Arguments

- formula:

  formule façon `lm`.

- data:

  data.frame.

- tau:

  niveau de quantile (défaut 0.5 = médiane / LAD).

- maxit:

  itérations maximales.

- tol:

  tolérance d'arrêt (variation des coefficients).

- eps:

  perturbation de la majoration (évite la division par 0).

## Value

liste : `coefficients`, `tau`, `fitted`, `residuals`, `loss`, `iter`.
