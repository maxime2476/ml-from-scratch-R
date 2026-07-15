# Batch normalization (passe avant)

Normalise chaque colonne (feature) sur le batch : \\\hat
x=(x-\mu)/\sqrt{ \sigma^2+\epsilon}\\, puis \\y=\gamma\hat x+\beta\\.
Stabilise et accelere l'apprentissage (reduit le decalage de covariance
interne).

## Usage

``` r
batch_norm(X, gamma = rep(1, ncol(X)), beta = rep(0, ncol(X)), eps = 1e-05)
```

## Arguments

- X:

  matrice n x d (batch).

- gamma, beta:

  parametres d'echelle et de decalage (longueur d).

- eps:

  stabilisateur.

## Value

liste : `out`, `cache` (pour la retropropagation).
