# Estimation ajustée pour un confondeur de force donnée (OVB)

Décale l'effet et l'erreur standard selon les R² partiels
**hypothétiques** du confondeur avec le traitement (`r2dz`) et avec le
résultat (`r2yz`).

## Usage

``` r
adjusted_estimate(estimate, se, df, r2dz, r2yz, reduce = TRUE)
```

## Arguments

- estimate:

  estimation de l'effet.

- se:

  erreur standard.

- df:

  degrés de liberté résiduels.

- r2dz:

  R² partiel confondeur-traitement.

- r2yz:

  R² partiel confondeur-résultat.

- reduce:

  réduire l'effet vers 0 (défaut TRUE).

## Value

liste : `estimate`, `se`, `t`, `bias`.
