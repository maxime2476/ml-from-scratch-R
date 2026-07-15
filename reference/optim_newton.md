# Newton-Raphson (minimisation)

Itère \\x\_{k+1} = x_k - \[\nabla^2 f(x_k)\]^{-1}\nabla f(x_k)\\ (éq.
0.17). Convergence quadratique locale (Th. 0.7). Le pas de Newton est
obtenu en résolvant \\\nabla^2 f\\\delta = \nabla f\\.

## Usage

``` r
optim_newton(grad, hess, x0, max_iter = 100L, tol = 1e-10, f = NULL)
```

## Arguments

- grad:

  fonction gradient `grad(x)`.

- hess:

  fonction hessienne `hess(x)` (matrice SPD localement).

- x0:

  point initial.

- max_iter:

  nombre maximal d'itérations.

- tol:

  seuil d'arrêt sur \\\\x\_{k+1}-x_k\\\\.

- f:

  (optionnel) fonction objectif.

## Value

liste : `par`, `iter`, `grad_norm`, `value`.
