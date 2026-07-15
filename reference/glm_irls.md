# Ajustement d'un GLM par IRLS

Boucle IRLS (éq. 3.7-3.9) : à chaque étape, poids \\W\\, réponse de
travail \\z\\ (éq. 3.8), puis WLS de z sur X (QR du Module 0).
Convergence sur la variation relative de la déviance (critère de `glm`).
Variance par l'information de Fisher \\(X^TWX)^{-1}\\ (éq. 3.10).

## Usage

``` r
glm_irls(
  formula,
  data,
  family = c("binomial", "poisson"),
  maxit = 25L,
  epsilon = 1e-08
)
```

## Arguments

- formula:

  formule façon `glm`.

- data:

  data.frame.

- family:

  "binomial" (logit) ou "poisson" (log).

- maxit:

  itérations maximales (défaut 25, comme `glm`).

- epsilon:

  tolérance de convergence (défaut 1e-8, comme `glm`).

## Value

objet `glm_irls` : `coefficients`, `vcov`, `se`, `fitted`,
`linear.predictors`, `deviance`, `null.deviance`, `loglik`, `iter`,
`df.residual`, `rank`, `family`, `weights`, `model_matrix`, `response`.
