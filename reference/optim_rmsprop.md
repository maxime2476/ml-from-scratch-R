# Optimiseur RMSprop

\\v_t=\rho v\_{t-1}+(1-\rho)g_t^2\\, \\x\_{t+1}=x_t-\eta
g_t/(\sqrt{v_t}+\epsilon)\\.

## Usage

``` r
optim_rmsprop(
  grad,
  x0,
  lr = 0.01,
  rho = 0.9,
  eps = 1e-08,
  max_iter = 5000L,
  tol = 1e-08
)
```

## Arguments

- grad, x0, lr:

  cf. `optim_adam` ; @param rho decroissance ; @param eps,max_iter,tol
  arret.

## Value

liste : `par`, `iters`.
