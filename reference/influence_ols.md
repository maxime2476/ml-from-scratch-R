# Fonction d'influence de l'OLS

\\\mathrm{IC}\_i = n\\(X^\top X)^{-1} x_i \hat\varepsilon_i\\, de sorte
que \\\hat\beta-\beta \approx \frac1n\sum_i \mathrm{IC}\_i\\ et
\\\widehat{\operatorname{Var}}(\hat\beta)=\frac1{n^2}\sum_i
\mathrm{IC}\_i\mathrm{IC}\_i^\top\\ — **exactement** le sandwich HC0
(Module 2).

## Usage

``` r
influence_ols(X, y)
```

## Arguments

- X:

  matrice de design n x p (constante incluse).

- y:

  réponse.

## Value

liste : `ic` (n x p), `vcov` (= sandwich HC0), `beta`.
