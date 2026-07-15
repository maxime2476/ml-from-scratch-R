# Retropropagation du max-pooling

Route le gradient vers la position du maximum de chaque bloc (les autres
recoivent 0).

## Usage

``` r
max_pool2d_backward(dout, cache)
```

## Arguments

- dout:

  gradient en sortie (oh x ow) ; @param cache sortie de `max_pool2d`.

## Value

`dX` (matrice de la taille de l'entree).
