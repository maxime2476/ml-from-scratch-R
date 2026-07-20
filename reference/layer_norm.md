# Normalisation de couche (layer norm)

Normalise chaque LIGNE (par exemple) puis remet a l'echelle : brique des
blocs de Transformer (avec les connexions residuelles).

## Usage

``` r
layer_norm(X, gamma = rep(1, ncol(X)), beta = rep(0, ncol(X)), eps = 1e-05)
```

## Arguments

- X:

  matrice

- gamma, beta:

  echelle et decalage (longueur ncol)

- eps:

  stabilisateur.

## Value

matrice normalisee.
