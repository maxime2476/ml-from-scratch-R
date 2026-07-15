# Prédiction conforme par découpage (split conformal)

Ajuste un modèle sur l'échantillon d'entraînement, calibre les scores de
non-conformité \\\|y-\hat\mu(x)\|\\ sur l'échantillon de calibration, et
renvoie des intervalles \\\hat\mu(x)\pm\hat q\\ pour `X_test` (Théorème
19.1). Avec `normalize = TRUE`, les scores sont divisés par une
dispersion locale \\\hat\sigma(x)\\ (intervalles de largeur variable).

## Usage

``` r
conformal_split(
  X_train,
  y_train,
  X_cal,
  y_cal,
  X_test,
  fit_fn,
  predict_fn,
  alpha = 0.1,
  normalize = FALSE,
  sigma_fn = NULL
)
```

## Arguments

- X_train, y_train:

  données d'entraînement.

- X_cal, y_cal:

  données de calibration.

- X_test:

  points où prédire.

- fit_fn:

  fonction `(X, y) -> modèle`.

- predict_fn:

  fonction `(modèle, X) -> prédictions`.

- alpha:

  niveau (défaut 0.1 -\> 90 %).

- normalize:

  normalisation locale (défaut FALSE).

- sigma_fn:

  (si normalize) fonction `(X) -> dispersion locale positive`.

## Value

liste : `lower`, `upper`, `pred`, `qhat`.
