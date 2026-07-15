# Rétropropagation : gradients analytiques (éq. 12.2-12.5)

Rétropropagation : gradients analytiques (éq. 12.2-12.5)

## Usage

``` r
mlp_backward(params, X, Y, activation, loss)
```

## Arguments

- params:

  liste des paramètres.

- X:

  matrice n x d0.

- Y:

  cible n x d2.

- activation:

  activation cachée.

- loss:

  "mse" ou "logloss".

## Value

liste des gradients `W1, b1, W2, b2`.
