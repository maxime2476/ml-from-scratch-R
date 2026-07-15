# Concentration des distances (éq. 7.4-7.5)

Pour n points i.i.d. en dimension p, calcule, pour chaque point de
requête, le rapport de contraste \\(D\_{\max}-D\_{\min})/D\_{\min}\\
entre distance la plus grande et la plus petite aux autres points, et le
coefficient de variation des distances au carré (théorique \\\propto
1/\sqrt p\\).

## Usage

``` r
distance_concentration(
  n,
  p,
  gen = function(n, p) matrix(runif(n * p), n, p),
  seed = NULL
)
```

## Arguments

- n:

  nombre de points.

- p:

  dimension.

- gen:

  générateur de points (défaut : uniforme sur `[0,1]^p`).

- seed:

  graine.

## Value

liste : `contrast` (moyenne de (Dmax-Dmin)/Dmin), `cv_d2` (coef. de
variation des distances au carré).
