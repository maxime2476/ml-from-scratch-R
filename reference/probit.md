# Probit (reponse binaire) par IRLS (lien probit)

\\P(y=1)=\Phi(X\beta)\\. Maximum de vraisemblance par moindres carres
ponderes iteres (IRLS, Module 3) avec lien probit : poids
\\w=\phi^2/\[\Phi(1-\Phi)\]\\, reponse de travail
\\z=\eta+(y-\Phi)/\phi\\.

## Usage

``` r
probit(formula, data, maxit = 50L)
```

## Arguments

- formula:

  formule

- data:

  data.frame

- maxit:

  iterations max

## Value

liste : `coefficients`, `vcov`, `se`, `fitted`, `loglik`.
