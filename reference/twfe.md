# Estimateur TWFE (two-way fixed effects) de l'effet du traitement

Régression de \\y\\ sur l'indicatrice \\D\_{it}\\ avec effets fixes
unité et temps (double démean, FWL — Modules 1, 21). SE groupée par
unité.

## Usage

``` r
twfe(data, yname = "y", idname = "id", tname = "t", gname = "g")
```

## Arguments

- data:

  data.frame du panel.

- yname, idname, tname, gname:

  noms des colonnes (résultat, unité, temps, cohorte de premier
  traitement ; jamais-traités = Inf ou 0).

## Value

liste : `coef`, `se` (groupée par unité).
