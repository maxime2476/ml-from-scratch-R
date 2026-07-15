# Tests t robustes (coeftest) avec une variance donnée

Reproduit
[`lmtest::coeftest`](https://rdrr.io/pkg/lmtest/man/coeftest.html) : t =
coef / se_robuste, loi de Student à `df.residual` degrés de liberté.

## Usage

``` r
coeftest_hc(fit, vcov = vcov_hc(fit, "HC3"))
```

## Arguments

- fit:

  objet `ols`.

- vcov:

  matrice de variance (p.ex. sortie de `vcov_hc`).

## Value

data.frame estimate/se/t/p_value.
