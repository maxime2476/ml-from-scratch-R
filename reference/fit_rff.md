# Ajuste un modèle à caractéristiques aléatoires (interpolant ou régularisé)

`lambda = 0` : moindres carrés de **norme minimale** (`solve_ls_svd`,
M0) — interpole les données dès que \\D \ge n\\. `lambda > 0` : ridge
(M4) sur les caractéristiques. Renvoie un prédicteur.

## Usage

``` r
fit_rff(X, y, D, gamma = 1, seed = 1, lambda = 0)
```

## Arguments

- X, y:

  données d'apprentissage.

- D:

  dimension des caractéristiques.

- gamma, seed:

  cf. `random_features`.

- lambda:

  pénalité ridge (0 = interpolation de norme minimale).

## Value

fonction `newX -> prédictions`.
