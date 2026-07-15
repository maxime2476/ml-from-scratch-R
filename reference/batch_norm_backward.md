# Batch normalization (passe arriere)

Gradients de la perte par rapport a \\X\\, \\\gamma\\, \\\beta\\ etant
donne `dout` (gradient en sortie).

## Usage

``` r
batch_norm_backward(dout, cache)
```

## Arguments

- dout:

  gradient en sortie (n x d).

- cache:

  sortie de `batch_norm`.

## Value

liste : `dX`, `dgamma`, `dbeta`.
