# Passe arrière : calcule tous les gradients par une seule rétropropagation

Amorce \\\partial L/\partial L=1\\ sur le nœud de sortie, puis parcourt
la bande en ordre inverse (topologique) en propageant chaque gradient à
ses parents. Après appel, `node$grad` de chaque `adnode` contient
\\\partial L/ \partial\text{node}\\.

## Usage

``` r
backward(node)
```

## Arguments

- node:

  `adnode` scalaire de sortie (la « perte »).
