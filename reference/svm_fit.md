# SVM a marge souple par le probleme dual (programmation quadratique)

Maximise \\\sum_i\alpha_i-\tfrac12\sum\_{ij}\alpha_i\alpha_j y_iy_j
K(x_i,x_j)\\ sous \\0\le\alpha_i\le C\\ et \\\sum_i\alpha_i y_i=0\\.
Resolu par
[`quadprog::solve.QP`](https://rdrr.io/pkg/quadprog/man/solve.QP.html).
Seuls les **vecteurs de support** (\\\alpha_i\>0\\) comptent.

## Usage

``` r
svm_fit(X, y, C = 1, kernel = c("linear", "rbf"), gamma = 0.5)
```

## Arguments

- X:

  matrice n x p

- y:

  etiquettes dans \\\\-1,+1\\\\.

- C:

  penalite (marge souple)

- kernel:

  "linear" ou "rbf"

- gamma:

  parametre RBF.

## Value

objet `svm` : `alpha`, `b`, vecteurs de support, hyperparametres.
