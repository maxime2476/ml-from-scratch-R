# Entraînement d'un MLP par SGD (éq. 12.6)

Entraînement d'un MLP par SGD (éq. 12.6)

## Usage

``` r
mlp_fit(
  X,
  y,
  hidden = 8L,
  activation = "tanh",
  loss = c("mse", "logloss"),
  epochs = 200L,
  lr = 0.05,
  batch = 32L,
  seed = NULL
)
```

## Arguments

- X:

  matrice n x d0.

- y:

  cible (vecteur ou matrice n x d2).

- hidden:

  largeur de la couche cachée d1.

- activation:

  "tanh" (défaut), "relu" ou "sigmoid".

- loss:

  "mse" ou "logloss".

- epochs:

  nombre d'époques.

- lr:

  taux d'apprentissage.

- batch:

  taille de mini-lot.

- seed:

  graine (initialisation + mélange).

## Value

objet `mlp` : `params`, `activation`, `loss`, `loss_hist`, `d`.
