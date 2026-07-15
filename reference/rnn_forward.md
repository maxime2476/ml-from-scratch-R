# Passe avant d'un RNN simple (sequence -\> sequence)

\\h_t=\tanh(W\_{xh}x_t+W\_{hh}h\_{t-1}+b_h)\\, \\y_t=W\_{hy}h_t+b_y\\.

## Usage

``` r
rnn_forward(X, Wxh, Whh, Why, bh, by, h0 = rep(0, nrow(Wxh)))
```

## Arguments

- X:

  sequence d'entree (T x d).

- Wxh, Whh, Why:

  matrices de poids (H x d), (H x H), (O x H).

- bh, by:

  biais (H), (O).

- h0:

  etat initial (defaut 0).

## Value

liste : `Y` (T x O), `H` (T x H), `cache`.
