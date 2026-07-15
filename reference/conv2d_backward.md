# Retropropagation de la convolution 2D

Gradients par rapport a l'entree, aux noyaux et aux biais, etant donne
`dout`.

## Usage

``` r
conv2d_backward(dout, cache)
```

## Arguments

- dout:

  gradient en sortie (Hout x Wout x F).

- cache:

  sortie de `conv2d`.

## Value

liste : `dX`, `dK`, `db`.
