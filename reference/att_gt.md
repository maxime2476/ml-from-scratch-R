# ATT groupe-temps (Callaway & Sant'Anna 2021)

\\\mathrm{ATT}(g,t)=\bigl\[\bar Y\_{g,t}-\bar Y\_{g,g-1}\bigr\]-
\bigl\[\bar Y\_{C,t}-\bar Y\_{C,g-1}\bigr\]\\, base universelle \\g-1\\,
contrôles = jamais-traités (`"never"`) ou pas-encore-traités
(`"notyet"`). Évite les « comparaisons interdites » du TWFE.

## Usage

``` r
att_gt(
  data,
  yname = "y",
  idname = "id",
  tname = "t",
  gname = "g",
  control = c("never", "notyet")
)
```

## Arguments

- data:

  data.frame du panel ; noms comme `twfe`.

- yname, idname, tname, gname:

  noms des colonnes (résultat, unité, temps, cohorte de premier
  traitement ; jamais-traités = Inf ou 0).

- control:

  "never" (jamais-traités) ou "notyet" (pas encore traités en t).

## Value

data.frame : `g`, `t`, `att`, `n_treat`.
