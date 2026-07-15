# Noyau gaussien (RBF)

\\k(x,x')=\sigma_f^2\exp(-\\x-x'\\^2/(2\ell^2))\\.

## Usage

``` r
rbf_kernel(X1, X2, lengthscale = 1, variance = 1)
```

## Arguments

- X1, X2:

  matrices (n1 x p), (n2 x p).

- lengthscale:

  échelle \\\ell\\.

- variance:

  variance du signal \\\sigma_f^2\\.

## Value

matrice de noyau n1 x n2.
