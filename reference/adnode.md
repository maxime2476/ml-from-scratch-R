# Crée un nœud de calcul (variable enregistrée sur la bande)

Un `adnode` est un environnement (sémantique de référence) portant sa
`value`, son `grad` accumulé (même forme), et une fonction `backward`
qui propage le gradient vers ses parents (règle de la chaîne).

## Usage

``` r
adnode(value)
```

## Arguments

- value:

  scalaire, vecteur ou matrice.

## Value

objet `adnode`.
