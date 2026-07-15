# Optimiseur Adam (moments adaptatifs)

Combine momentum (moment 1) et mise a l'echelle par la variance des
gradients (moment 2), avec correction de biais : \\m_t=\beta_1
m\_{t-1}+(1-\beta_1)g_t\\, \\v_t=\beta_2 v\_{t-1}+(1-\beta_2)g_t^2\\,
pas \\x\_{t+1}=x_t-\eta\\\hat m_t/(\sqrt{\hat v_t}+\epsilon)\\.

## Usage

``` r
optim_adam(
  grad,
  x0,
  lr = 0.01,
  beta1 = 0.9,
  beta2 = 0.999,
  eps = 1e-08,
  max_iter = 5000L,
  tol = 1e-08
)
```

## Arguments

- grad:

  fonction `x -> gradient`.

- x0:

  point initial.

- lr:

  pas d'apprentissage.

- beta1, beta2:

  taux de decroissance des moments.

- eps:

  stabilisateur numerique.

- max_iter, tol:

  arret.

## Value

liste : `par`, `iters`, `path` (valeurs successives).
