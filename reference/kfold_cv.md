# Validation croisée K-fold générique (éq. 6.2)

Validation croisée K-fold générique (éq. 6.2)

## Usage

``` r
kfold_cv(
  X,
  y,
  K = 10L,
  fit_fun = function(Xtr, ytr) solve_ls_qr(Xtr, ytr)$coefficients,
  pred_fun = function(beta, Xte) as.numeric(Xte %*% beta),
  seed = NULL
)
```

## Arguments

- X:

  matrice de design (constante incluse si le modèle en a une).

- y:

  réponse.

- K:

  nombre de blocs (défaut 10).

- fit_fun:

  fonction `(Xtr, ytr) -> modèle` (défaut : OLS via QR).

- pred_fun:

  fonction `(modèle, Xte) -> prédictions` (défaut : OLS).

- seed:

  graine pour le tirage des blocs.

## Value

liste : `cv`, `se`, `fold_errors`, `folds`.
