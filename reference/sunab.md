# Event-study de Sun & Abraham (2021) — interaction-weighted

Régression saturée \\y\\ sur les indicatrices cohorte×ancienneté
(référence \\e=-1\\) avec effets fixes unité+temps, puis agrégation des
coefficients par ancienneté \\e\\ pondérée par les parts de cohorte.
Robuste aux effets dynamiques hétérogènes, contrairement au TWFE.

## Usage

``` r
sunab(data, yname = "y", idname = "id", tname = "t", gname = "g")
```

## Arguments

- data:

  data.frame du panel ; noms comme `twfe`.

- yname, idname, tname, gname:

  noms des colonnes (résultat, unité, temps, cohorte de premier
  traitement ; jamais-traités = Inf ou 0).

## Value

liste : `es` (data.frame `e`, `att` — event-study) et `att` (ATT global
post-traitement, pondéré par les tailles de cellule).
