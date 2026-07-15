# Taille d'echantillon effective (ESS)

\\\mathrm{ESS}=n/(1+2\sum_k\hat\rho_k)\\, somme des autocorrelations
tronquee au premier terme negatif (sequence positive initiale de Geyer).
Mesure combien de tirages INDEPENDANTS equivalent a la chaine correlee.

## Usage

``` r
ess(x)
```

## Arguments

- x:

  chaine (vecteur).

## Value

la taille d'echantillon effective.
