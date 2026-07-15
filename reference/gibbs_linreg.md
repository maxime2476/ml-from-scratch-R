# Gibbs pour la regression lineaire bayesienne (prior non informatif)

Alterne les tirages des lois CONDITIONNELLES conjuguees :
\\\beta\mid\sigma^2,y\sim\mathcal N(\hat\beta\_{OLS},\sigma^2(X^\top
X)^{-1})\\ et
\\\sigma^2\mid\beta,y\sim\text{Inv-Gamma}(n/2,\\\\y-X\beta\\^2/2)\\. La
loi stationnaire est le posterieur joint ; sa moyenne coincide avec
l'OLS (M1).

## Usage

``` r
gibbs_linreg(X, y, n_iter = 5000L, burn = 1000L)
```

## Arguments

- X:

  design

- y:

  reponse

- n_iter:

  iterations

- burn:

  rodage

## Value

liste : `beta` (echantillons apres burn), `sigma2`.
