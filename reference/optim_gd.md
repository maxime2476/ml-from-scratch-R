# Descente de gradient à pas constant

Itère \\x\_{k+1} = x_k - t\\\nabla f(x_k)\\ (éq. 0.15). Sous \\f\\
convexe L-lisse et pas \\t = 1/L\\, vitesse \\O(1/k)\\ (Th. 0.6).

## Usage

``` r
optim_gd(grad, x0, step, max_iter = 10000L, tol = 1e-08, f = NULL)
```

## Arguments

- grad:

  fonction gradient : `grad(x)` renvoie \\\nabla f(x)\\.

- x0:

  point initial.

- step:

  pas constant \\t\\ (idéalement \\1/L\\).

- max_iter:

  nombre maximal d'itérations.

- tol:

  seuil d'arrêt sur \\\\x\_{k+1}-x_k\\\\.

- f:

  (optionnel) fonction objectif, pour renvoyer la valeur finale.

## Value

liste : `par`, `iter`, `grad_norm`, `value`.
