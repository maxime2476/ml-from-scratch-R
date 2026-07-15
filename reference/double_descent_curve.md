# Courbe de risque en fonction de la complexité (double descente)

Balaie la dimension \\D\\ des caractéristiques et renvoie l'erreur
quadratique d'**apprentissage** et de **test** pour chaque \\D\\. Avec
`lambda = 0`, le test présente un **pic à \\D\approx n\\** (seuil
d'interpolation) puis **redescend** — la double descente.

## Usage

``` r
double_descent_curve(Xtr, ytr, Xte, yte, Ds, gamma = 1, seed = 1, lambda = 0)
```

## Arguments

- Xtr, ytr:

  apprentissage.

- Xte, yte:

  test.

- Ds:

  vecteur des dimensions à essayer.

- gamma, seed:

  cf. `random_features`.

- lambda:

  pénalité ridge.

## Value

data.frame : `D`, `train_mse`, `test_mse`.
