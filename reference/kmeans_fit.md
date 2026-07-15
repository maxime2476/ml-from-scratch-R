# k-means par l'algorithme de Lloyd (éq. 11.4)

Minimisation alternée de l'inertie intra-classe (affectation au centre
le plus proche, puis moyenne de classe). Reproduit
`stats::kmeans(algorithm="Lloyd")` à initialisation identique.

## Usage

``` r
kmeans_fit(X, K, centers = NULL, max_iter = 100L, nstart = 1L, seed = NULL)
```

## Arguments

- X:

  matrice n x p.

- K:

  nombre de classes.

- centers:

  matrice K x p de centres initiaux (sinon tirage aléatoire).

- max_iter:

  itérations maximales.

- nstart:

  nombre de redémarrages aléatoires (si centers non fourni).

- seed:

  graine.

## Value

liste : `cluster`, `centers`, `withinss`, `tot_withinss`, `iter`.
