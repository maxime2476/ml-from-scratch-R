# Encodage positionnel sinusoidal

\\PE\_{pos,2i}=\sin(pos/10000^{2i/d})\\, \\PE\_{pos,2i+1}=\cos(\cdot)\\.
Ajoute l'information d'ORDRE (que l'attention, permutation-equivariante,
ignore).

## Usage

``` r
positional_encoding(seq_len, d_model)
```

## Arguments

- seq_len:

  longueur de la sequence

- d_model:

  dimension du modele.

## Value

matrice (seq_len x d_model).
