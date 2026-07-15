# Regression locale lineaire

Ajuste en chaque \\x_0\\ une droite ponderee par le noyau (WLS locale)
et renvoie l'ordonnee a l'origine. Corrige le **biais de bord** du
Nadaraya-Watson (biais d'ordre \\h^2\\ uniforme jusqu'aux bords).

## Usage

``` r
local_linear(x, y, x0, bw, kernel = c("gauss", "tri"))
```

## Arguments

- x, y:

  donnees

- x0:

  points d'evaluation

- bw:

  fenetre

- kernel:

  "gauss" ou "tri" (triangulaire).

## Value

vecteur des valeurs ajustees en `x0`.
