# Analyse discriminante lineaire (LDA)

Modele : \\x\mid y=k\sim\mathcal N(\mu_k,\Sigma)\\ (covariance
**commune**). La regle de Bayes donne une frontiere **lineaire** ; le
discriminant est
\\\delta_k(x)=x^\top\Sigma^{-1}\mu_k-\tfrac12\mu_k^\top\Sigma^{-1}\mu_k+\log\pi_k\\.

## Usage

``` r
lda_fit(X, y)
```

## Arguments

- X:

  matrice n x p ; @param y etiquettes de classe.

## Value

objet `lda_scratch`.
