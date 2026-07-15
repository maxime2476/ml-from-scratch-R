# Estimateur « en un pas » (one-step / Newton-scoring)

\\\tilde\theta = \hat\theta_0 + \frac1n\sum_i
\mathrm{IC}\_i(\hat\theta_0)\\ : un estimateur initial grossier corrigé
de sa fonction d'influence devient **asymptotiquement efficace**. C'est
la forme abstraite de la correction du lasso débiaisé (Module 22) et du
score orthogonal de Neyman (Module 16).

## Usage

``` r
onestep(theta0, ic_at)
```

## Arguments

- theta0:

  estimateur initial (scalaire ou vecteur).

- ic_at:

  fonction `theta -> IC(theta)` : vecteur de longueur n (theta scalaire)
  ou matrice n x p (theta vectoriel).

## Value

l'estimateur corrigé \\\tilde\theta\\.
