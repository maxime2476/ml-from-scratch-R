# Classifieur naif de Bayes (gaussien)

Hypothese **naive** : les features sont independantes conditionnellement
a la classe, \\P(x\mid y=k)=\prod_j\mathcal
N(x_j;\mu\_{kj},\sigma\_{kj}^2)\\. Rapide, robuste en haute dimension,
mais biaise si les features sont correlees.

## Usage

``` r
naive_bayes_fit(X, y)
```

## Arguments

- X, y:

  cf. `lda_fit`.

## Value

objet `nb_scratch`.
