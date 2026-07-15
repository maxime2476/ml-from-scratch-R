# Estimation d'un ARMA(p,q) par moindres carres conditionnels (CSS)

Minimise \\\sum_t \hat\varepsilon_t^2\\, ou
\\\hat\varepsilon_t=(x_t-\mu)- \sum
\phi_i(x\_{t-i}-\mu)-\sum\theta_j\hat\varepsilon\_{t-j}\\ est calcule
recursivement (erreurs initiales nulles). Optimisation via `optim`.

## Usage

``` r
arma_css(x, p = 1L, q = 1L)
```

## Arguments

- x:

  serie

- p, q:

  ordres AR et MA

## Value

liste : `ar`, `ma`, `mean`, `sigma2`.
