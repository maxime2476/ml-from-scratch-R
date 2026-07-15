# Opérateur de soft-thresholding (éq. 4.9)

\\\mathcal S(z,\lambda) = \mathrm{sign}(z)\\(\|z\|-\lambda)\_+\\.
Vectorisé.

## Usage

``` r
soft_threshold(z, lambda)
```

## Arguments

- z:

  scalaire ou vecteur.

- lambda:

  seuil (\>= 0).

## Value

la (les) valeur(s) seuillée(s).
