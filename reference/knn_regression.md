# Régression KNN (éq. 7.1)

Prédit en chaque point de test la moyenne des réponses de ses k plus
proches voisins d'apprentissage.

## Usage

``` r
knn_regression(X_train, y_train, X_test, k)
```

## Arguments

- X_train:

  matrice n x p d'apprentissage.

- y_train:

  réponses d'apprentissage (numériques).

- X_test:

  matrice m x p de test.

- k:

  nombre de voisins.

## Value

vecteur des m prédictions.
