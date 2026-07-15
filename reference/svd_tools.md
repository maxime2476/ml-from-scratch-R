# Outils SVD : rang numérique, conditionnement, pseudo-inverse

À partir de \\X = U\Sigma V^T\\ (éq. 0.11), calcule le rang numérique,
le conditionnement \\\kappa_2 = \sigma\_{\max}/\sigma\_{\min}\\ (éq.
0.3) et la pseudo-inverse de Moore-Penrose \\X^+ = V\Sigma^+ U^T\\ (éq.
0.12). La SVD elle-même est déléguée à
[`svd()`](https://rdrr.io/r/base/svd.html) (existence admise, Th. 0.4).

## Usage

``` r
svd_tools(X, tol = NULL)
```

## Arguments

- X:

  matrice n x p.

- tol:

  seuil de rang ; par défaut \\\max(n,p)\\\varepsilon\\\sigma_1\\.

## Value

liste : `d` (valeurs singulières), `rank`, `kappa`, `pinv`
(pseudo-inverse p x n), `u`, `v`, `tol`.
