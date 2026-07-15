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

  matrice n x p ; @param dims dimension de sortie (2) ; @param
  perplexity ;

- iter, eta, sigma:

  nb d'iterations, pas, largeur des gaussiennes.

## Value

matrice n x dims (le plongement).
