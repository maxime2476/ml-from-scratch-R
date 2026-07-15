# Descente de gradient a momentum (boule pesante de Polyak)

\\u_t=\mu u\_{t-1}-\eta g_t\\, \\x\_{t+1}=x_t+u_t\\.

## Usage

``` r
optim_momentum(
  grad,
  x0,
  step = 0.01,
  momentum = 0.9,
  max_iter = 5000L,
  tol = 1e-08
)
```

## Arguments

- grad, x0:

  cf. `optim_adam` ; @param step pas ; @param momentum \\\mu\\ ;

- max_iter, tol:

  arret.

## Value

liste : `par`, `iters`.
