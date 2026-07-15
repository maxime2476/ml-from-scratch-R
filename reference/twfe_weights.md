# Poids de la régression TWFE (de Chaisemartin & D'Haultfœuille 2020)

Le coefficient TWFE est \\\sum\_{(i,t):D=1} w\_{it}\\\tau\_{it}\\ avec
\\w\_{it}=\tilde D\_{it}/\sum \tilde D^2\\ (\\\tilde D\\ = traitement
résidualisé des effets fixes). Les poids somment à 1 mais certains sont
**négatifs** : TWFE n'est pas une moyenne convexe des effets — il peut
même être de signe opposé à tous les effets individuels.

## Usage

``` r
twfe_weights(data, yname = "y", idname = "id", tname = "t", gname = "g")
```

## Arguments

- data:

  data.frame du panel ; noms de colonnes comme `twfe`.

- yname, idname, tname, gname:

  noms des colonnes (résultat, unité, temps, cohorte de premier
  traitement ; jamais-traités = Inf ou 0).

## Value

liste : `weights` (sur cellules traitées), `share_negative`, `sum`.
