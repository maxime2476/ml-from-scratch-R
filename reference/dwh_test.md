# Test d'endogeneite de Durbin-Wu-Hausman (regression augmentee)

Regresse les regresseurs suspects sur les instruments, recupere les
residus, les ajoute a l'OLS et teste leur significativite jointe (F).
Rejet = OLS biaise, l'IV (Module 5) est requis.

## Usage

``` r
dwh_test(y, X, Z, endog)
```

## Arguments

- y:

  reponse

- X:

  design (constante + regresseurs, endogenes inclus)

- Z:

  instruments (constante + exogenes + exclus).

- endog:

  indices des colonnes ENDOGENES de X.

## Value

liste : `statistic` (F), `df1`, `df2`, `p_value`.
