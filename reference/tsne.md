# t-SNE (visualisation preservant le voisinage) — version compacte

Convertit les distances en probabilites de voisinage (gaussiennes en
haute dimension, Student-t en basse) et minimise la divergence KL par
descente de gradient. Preserve la structure LOCALE : les clusters se
separent nettement.

## Usage

``` r
tsne(X, dims = 2L, perplexity = 30, iter = 500L, eta = 200, sigma = NULL)
```

## Arguments

- X:

  matrice n x p

- dims:

  dimension de sortie (2)

- perplexity:

  nombre effectif de voisins vise (eq. 42.5) ; fixe un sigma_i PAR POINT
  par dichotomie. Ignore si `sigma` est fourni.

- iter, eta:

  nb d'iterations, pas de la descente.

- sigma:

  largeur commune imposee a tous les points ; par defaut `NULL`, et les
  sigma_i sont calibres depuis `perplexity`.

## Value

matrice n x dims (le plongement).
