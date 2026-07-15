# Retropropagation dans le temps (BPTT) du RNN simple

Retropropagation dans le temps (BPTT) du RNN simple

## Usage

``` r
rnn_backward(dY, cache)
```

## Arguments

- dY:

  gradient en sortie (T x O).

- cache:

  sortie de `rnn_forward`.

## Value

liste : `dWxh`, `dWhh`, `dWhy`, `dbh`, `dby`.
