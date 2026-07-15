# Variance sandwich d'un M-estimateur (éq. 14.2)

\\\hat V = \hat A^{-1}\hat B\hat A^{-1}/n\\, avec \\\hat B =
\frac1n\sum_i \psi_i\psi_i^\top\\ et \\\hat A\\ fournie (moyenne des
dérivées de la fonction d'estimation, \\A = \mathbb
E\[-\partial\_\theta\psi\]\\).

## Usage

``` r
m_estimation_vcov(psi, A)
```

## Arguments

- psi:

  matrice n x p des valeurs de la fonction d'estimation \\\psi_i\\.

- A:

  matrice p x p, \\\hat A\\ (p.ex. \\X^\top X/n\\ pour l'OLS).

## Value

la matrice de variance sandwich p x p de \\\hat\theta\\.
