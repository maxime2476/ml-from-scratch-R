# Passe avant du MLP (éq. 12.1)

Passe avant du MLP (éq. 12.1)

## Usage

``` r
mlp_forward(params, X, activation)
```

## Arguments

- params:

  liste `W1, b1, W2, b2`.

- X:

  matrice n x d0.

- activation:

  activation de la couche cachée ("tanh", "relu", "sigmoid").

## Value

liste `Z1, A1, Z2`.
