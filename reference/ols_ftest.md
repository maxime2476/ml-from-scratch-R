# Test F de restrictions linéaires R beta = r (éq. 1.9)

Test F de restrictions linéaires R beta = r (éq. 1.9)

## Usage

``` r
ols_ftest(object, R, r = rep(0, nrow(R)))
```

## Arguments

- object:

  objet `ols`.

- R:

  matrice q x p des restrictions (rang q).

- r:

  vecteur cible (longueur q ; défaut zéro).

## Value

liste : `F`, `q`, `df.residual`, `p_value`.
