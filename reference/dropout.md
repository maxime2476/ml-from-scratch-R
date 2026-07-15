# Dropout (inverted dropout)

A l'apprentissage, met a zero une fraction `rate` des unites et remet a
l'echelle par \\1/(1-\text{rate})\\ (l'esperance est preservee). A
l'evaluation (`training=FALSE`), l'identite. Regularise en empechant la
co-adaptation des neurones.

## Usage

``` r
dropout(x, rate = 0.5, training = TRUE)
```

## Arguments

- x:

  activations (vecteur ou matrice).

- rate:

  probabilite d'extinction.

- training:

  TRUE (apprentissage) ou FALSE (evaluation).

## Value

liste : `out`, `mask`.
