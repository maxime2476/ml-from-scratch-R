# Lasso par coordinate descent (éq. 4.10)

Boucle de descente par coordonnées avec mise à jour du résidu ; chaque
coordonnée est mise à jour par soft-thresholding (Prop. 4.2).
L'intercept n'est pas pénalisé.

## Usage

``` r
lasso_fit(
  X,
  y,
  lambda,
  standardize = TRUE,
  intercept = TRUE,
  maxit = 10000L,
  tol = 1e-09
)
```

## Arguments

- X:

  matrice de design n x p (sans constante).

- y:

  réponse.

- lambda:

  pénalité \\\lambda \ge 0\\ (convention \\\tfrac12\\y-X\beta\\^2 +
  \lambda\\\beta\\\_1\\).

- standardize:

  centrer-réduire X (défaut TRUE).

- intercept:

  intercept non pénalisé (défaut TRUE).

- maxit:

  balayages maximaux.

- tol:

  tolérance d'arrêt (variation max des coefficients).

## Value

liste : `coefficients`, `beta`, `intercept`, `lambda`, `iter`, `fitted`.
