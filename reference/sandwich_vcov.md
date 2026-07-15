# Matrice de variance sandwich générique (éq. 2.1)

Calcule \\(X^TX)^{-1}(X^T \Omega X)(X^TX)^{-1}\\ pour une « viande »
\\X^T \Omega X\\ fournie. Sert de brique commune à HC et Newey-West.

## Usage

``` r
sandwich_vcov(fit, meat)
```

## Arguments

- fit:

  objet `ols`.

- meat:

  matrice p x p (la viande \\X^T \Omega X\\ estimée).

## Value

matrice de variance p x p.
