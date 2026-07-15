# Ajoute une colonne de 1 (terme de biais) — opération enregistrée

Pratique pour un MLP : \\\[\\A\\\mathbf 1\\\]\\. Le gradient est
restreint aux colonnes d'origine.

## Usage

``` r
ad_cbind1(a)
```

## Arguments

- a:

  `adnode` matrice n x k.

## Value

`adnode` matrice n x (k+1).
