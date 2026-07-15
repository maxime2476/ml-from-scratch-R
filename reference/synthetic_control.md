# Controle synthetique (Abadie)

Construit un "contrefactuel" de l'unite traitee comme combinaison
**convexe** des unites de controle (poids \\w\ge0\\, \\\sum w=1\\) qui
reproduit au mieux la trajectoire PRE-traitement. L'effet est l'ecart
post-traitement entre l'unite traitee et son synthetique. Poids par
programmation quadratique.

## Usage

``` r
synthetic_control(Y1, Y0, pre)
```

## Arguments

- Y1:

  trajectoire de l'unite traitee (longueur T).

- Y0:

  matrice T x J des unites de controle (donor pool).

- pre:

  indices des periodes PRE-traitement.

## Value

liste : `weights`, `synthetic` (trajectoire), `effect`
(post-traitement).
