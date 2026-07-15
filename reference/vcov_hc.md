# Variances robustes à l'hétéroscédasticité HC0–HC3 (éq. 2.2)

Viande \\\sum_i \phi_i \hat\varepsilon_i^2 x_i x_i^T\\ avec le facteur
\\\phi_i\\ dépendant du type (Prop. 2.2 pour la correction de levier).

## Usage

``` r
vcov_hc(fit, type = c("HC3", "HC0", "HC1", "HC2"))
```

## Arguments

- fit:

  objet `ols`.

- type:

  "HC0", "HC1", "HC2" ou "HC3".

## Value

matrice de variance p x p (compatible
[`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)).
