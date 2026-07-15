# Factorisation en matrices non negatives (NMF)

\\V\approx WH\\ avec \\W,H\ge 0\\ (mises a jour multiplicatives de
Lee-Seung). La contrainte de non-negativite donne une decomposition en
**parts additives** interpretables (contrairement a la PCA, qui peut
soustraire).

## Usage

``` r
nmf(V, k, iter = 500L)
```

## Arguments

- V:

  matrice non negative n x m ; @param k rang ; @param iter iterations.

## Value

liste : `W` (n x k), `H` (k x m), `reconstruction`.
