# L-BFGS (quasi-Newton à mémoire limitée)

Approxime l'inverse de la hessienne à partir des \\m\\ derniers couples
\\(s_k,y_k)\\ par la **récursion à deux boucles**, sans stocker de
matrice \\d\times d\\. Recherche linéaire d'Armijo (rétrogression) si
`f` est fourni. Convergence super-linéaire, coût par itération
\\O(md)\\.

## Usage

``` r
optim_lbfgs(grad, x0, f = NULL, m = 10L, max_iter = 200L, tol = 1e-08)
```

## Arguments

- grad:

  fonction gradient `grad(x)`.

- x0:

  point initial.

- f:

  fonction objectif (pour la recherche linéaire ; recommandé).

- m:

  taille de mémoire (défaut 10).

- max_iter:

  itérations maximales.

- tol:

  seuil sur la norme du gradient.

## Value

liste : `par`, `iter`, `grad_norm`, `value`.
