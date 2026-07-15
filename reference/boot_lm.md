# Bootstrap d'une régression linéaire (pairs ou résidus)

"pairs" rééchantillonne les couples \\(x_i,y_i)\\ (robuste à
l'hétéroscédasticité ; analogue du sandwich) ; "residual"
rééchantillonne les résidus (suppose des erreurs i.i.d.).

## Usage

``` r
boot_lm(formula, data, R = 2000L, method = c("pairs", "residual"), seed = NULL)
```

## Arguments

- formula:

  formule.

- data:

  data.frame.

- R:

  nombre de rééchantillons.

- method:

  "pairs" ou "residual".

- seed:

  graine.

## Value

objet `boot_lm` : `t0` (coefficients), `replicates` (R x p), `se`.
