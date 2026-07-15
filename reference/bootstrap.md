# Bootstrap non paramétrique d'une statistique

Rééchantillonne les données avec remise et applique `stat_fn` (principe
plug-in). Renvoie les répliques et l'erreur standard bootstrap (éq.
17.1).

## Usage

``` r
bootstrap(data, stat_fn, R = 2000L, seed = NULL)
```

## Arguments

- data:

  vecteur, matrice ou data.frame (les lignes sont rééchantillonnées).

- stat_fn:

  fonction `data -> statistique scalaire`.

- R:

  nombre de rééchantillons.

- seed:

  graine.

## Value

objet `bootstrap` : `t0`, `replicates`, `se`, `data`, `stat_fn`, `n`.
