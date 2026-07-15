# Regression Tobit (reponse censuree) par maximum de vraisemblance

Modele a variable latente \\y^\*=X\beta+\varepsilon\\, \\\varepsilon\sim
\mathcal N(0,\sigma^2)\\, observe \\y=\max(L,y^\*)\\. Vraisemblance :
densite normale pour les non censures, \\\Phi((L-X\beta)/\sigma)\\ pour
les censures.

## Usage

``` r
tobit_fit(formula, data, left = 0)
```

## Arguments

- formula:

  formule

- data:

  data.frame

- left:

  seuil de censure

## Value

liste : `coefficients`, `sigma`, `se`, `loglik`.
