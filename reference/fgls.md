# Moindres carres generalises FAISABLES (FGLS) pour heteroscedasticite

Modele de variance : \\\log\hat\varepsilon^2 = X\gamma\\. On en deduit
des poids \\\hat w_i = 1/\exp(X_i\hat\gamma)\\ et l'on refait une WLS
(Module 2). Plus efficace que l'OLS sous heteroscedasticite bien
specifiee.

## Usage

``` r
fgls(formula, data)
```

## Arguments

- formula, data:

  cf. `bp_test`.

## Value

liste : `coefficients`, `weights`, `se` (classiques WLS).
