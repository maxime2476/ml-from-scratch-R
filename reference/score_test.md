# Test du score / Rao (éq. 3.13) pour modèles emboîtés

Statistique de Rao = réduction de la somme des carrés pondérée obtenue
en régressant les résidus de travail du modèle réduit sur le design
complet, avec les poids de travail du réduit (une itération de score de
Fisher). C'est la forme équivalente à \\U(\tilde\beta)^T \mathcal
I(\tilde\beta)^{-1} U(\tilde\beta)\\ quand le réduit satisfait ses
équations normales, et c'est l'algorithme exact de
`anova.glm(test = "Rao")`.

## Usage

``` r
score_test(fit_full, fit_reduced)
```

## Arguments

- fit_full:

  objet `glm_irls` du modèle complet (fournit le design X).

- fit_reduced:

  objet `glm_irls` du modèle réduit (emboîté).

## Value

liste : `statistic`, `df`, `p_value`.
