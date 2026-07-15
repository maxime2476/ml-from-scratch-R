# Cp de Mallows (éq. 6.4)

\\C_p = \overline{\mathrm{err}} + 2\sigma^2 p/n\\, correction de
l'optimisme.

## Usage

``` r
mallows_cp(fit, sigma2)
```

## Arguments

- fit:

  objet `ols` (Module 1).

- sigma2:

  estimation de la variance du bruit (p.ex. d'un modèle riche).

## Value

la valeur du Cp (échelle erreur quadratique moyenne).
