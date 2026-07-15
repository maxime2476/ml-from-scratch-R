# Analyse discriminante quadratique (QDA)

Chaque classe a sa **propre** covariance \\\Sigma_k\\ : frontiere
**quadratique**.
\\\delta_k(x)=-\tfrac12\log\|\Sigma_k\|-\tfrac12(x-\mu_k)^\top
\Sigma_k^{-1}(x-\mu_k)+\log\pi_k\\.

## Usage

``` r
qda_fit(X, y)
```

## Arguments

- X, y:

  cf. `lda_fit`.

## Value

objet `qda_scratch`.
