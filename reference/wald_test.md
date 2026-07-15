# Test de Wald (éq. 3.11) pour H0 : R beta = r

Test de Wald (éq. 3.11) pour H0 : R beta = r

## Usage

``` r
wald_test(fit, R, r = rep(0, nrow(R)))
```

## Arguments

- fit:

  objet `glm_irls`.

- R:

  matrice q x p des restrictions.

- r:

  vecteur cible (défaut zéro).

## Value

liste : `statistic`, `df`, `p_value`.
