# Robustness value (Cinelli-Hazlett 2020)

Force minimale (en R² partiel) qu'un confondeur inobservé — supposé
aussi associé au traitement qu'au résultat — devrait avoir pour
**réduire l'effet de q×100 %** (`alpha = 1`), ou pour le rendre **non
significatif au seuil alpha** (`alpha < 1`).

## Usage

``` r
robustness_value(t, df, q = 1, alpha = 1)
```

## Arguments

- t:

  statistique t de l'effet.

- df:

  degrés de liberté résiduels.

- q:

  fraction de réduction visée (1 = ramener à 0).

- alpha:

  seuil de significativité (1 = point seulement).

## Value

la robustness value dans `[0,1]`.
