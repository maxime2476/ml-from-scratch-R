# Gradient accéléré de Nesterov

Ajoute un terme d'inertie (« momentum ») à la descente de gradient :
pour f convexe L-lisse et pas \\t=1/L\\, la vitesse passe de \\O(1/k)\\
(gradient) à \\O(1/k^2)\\ — l'accélération optimale au premier ordre
(Nesterov 1983).

## Usage

``` r
optim_nesterov(grad, x0, step, max_iter = 10000L, tol = 1e-08, f = NULL)
```

## Arguments

- grad:

  fonction gradient `grad(x)`.

- x0:

  point initial.

- step:

  pas \\t=1/L\\.

- max_iter:

  itérations maximales.

- tol:

  seuil d'arrêt sur \\\\x\_{k+1}-x_k\\\\.

- f:

  (optionnel) fonction objectif.

## Value

liste : `par`, `iter`, `grad_norm`, `value`.
