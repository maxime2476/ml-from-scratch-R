# Test du rapport de vraisemblance (éq. 3.12) pour modèles emboîtés

\\\mathrm{LR} = D\_{\text{réduit}} - D\_{\text{complet}}\\.

## Usage

``` r
lr_test(fit_full, fit_reduced)
```

## Arguments

- fit_full:

  objet `glm_irls` du modèle complet.

- fit_reduced:

  objet `glm_irls` du modèle réduit (emboîté).

## Value

liste : `statistic`, `df`, `p_value`.
