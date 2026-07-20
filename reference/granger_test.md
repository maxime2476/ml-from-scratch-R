# Test de causalite de Granger

\\x\\ **cause au sens de Granger** \\y\\ si ses retards ameliorent la
prediction de \\y\\. Test F comparant l'equation de \\y\\ avec et sans
les retards de \\x\\.

## Usage

``` r
granger_test(object, cause, effect)
```

## Arguments

- object:

  objet `var_fit`

- cause:

  indice de la variable causante

- effect:

  indice de la variable expliquee.

## Value

liste : `statistic` (F), `df1`, `df2`, `p_value`.
