# Fonction d'estimation de l'OLS : psi_i = x_i \* e_i

Fonction d'estimation de l'OLS : psi_i = x_i \* e_i

## Usage

``` r
ols_psi(X, resid)
```

## Arguments

- X:

  matrice de design n x p (constante incluse).

- resid:

  résidus \\\hat\varepsilon_i\\.

## Value

matrice n x p des \\\psi_i = x_i \hat\varepsilon_i\\.
