# Passe avant d'une cellule LSTM (sequence)

Portes d'entree/oubli/sortie et etat de cellule :
\\i=\sigma(W_i\[x_t,h\_{t-1}\])\\, \\f=\sigma(\cdot)\\,
\\o=\sigma(\cdot)\\, \\g=\tanh(\cdot)\\ ; \\c_t=f\odot c\_{t-1}+i\odot
g\\, \\h_t=o\odot\tanh(c_t)\\. L'etat de cellule \\c_t\\ forme une
"autoroute" ou le gradient circule sans s'evanouir.

## Usage

``` r
lstm_forward(X, Wi, Wf, Wo, Wg, bi, bf, bo, bg)
```

## Arguments

- X:

  sequence (T x d).

- Wi, Wf, Wo, Wg:

  matrices (H x (d+H)) empilant \\\[x_t,h\_{t-1}\]\\.

- bi, bf, bo, bg:

  biais (H).

## Value

liste : `H` (T x H), `C` (T x H).
