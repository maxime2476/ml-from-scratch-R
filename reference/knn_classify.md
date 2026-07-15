# Classification KNN par vote majoritaire

Vote majoritaire des k plus proches voisins. Les égalités de vote sont
rompues par la première classe (ordre de `levels`) — pour une
comparaison déterministe avec
[`class::knn`](https://rdrr.io/pkg/class/man/knn.html), utiliser k
impair et 2 classes (pas d'égalité).

## Usage

``` r
knn_classify(X_train, y_train, X_test, k)
```

## Arguments

- X_train:

  matrice n x p d'apprentissage.

- y_train:

  étiquettes d'apprentissage (facteur ou vecteur).

- X_test:

  matrice m x p de test.

- k:

  nombre de voisins.

## Value

facteur des m classes prédites (mêmes niveaux que y_train).
