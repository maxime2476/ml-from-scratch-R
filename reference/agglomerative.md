# Clustering hierarchique agglomeratif (assignation a k groupes)

Fusionne iterativement les deux groupes les plus proches selon la
**liaison** (mise a jour de Lance-Williams), jusqu'a `k` groupes.
Liaisons : "complete" (diametre), "single" (chaine), "average"
(moyenne).

## Usage

``` r
agglomerative(X, k, linkage = c("complete", "single", "average"))
```

## Arguments

- X:

  matrice n x p (ou objet `dist`)

- k:

  nombre de groupes

- linkage:

  "complete", "single" ou "average".

## Value

vecteur d'assignation (longueur n).
