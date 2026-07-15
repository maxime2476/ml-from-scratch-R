# Analyse en composantes principales par SVD (éq. 11.2-11.3)

Centre X puis calcule l'ACP via la SVD (voie 2 = voie 1 par Prop. 11.2).
Reproduit `prcomp` : `sdev` = d/sqrt(n-1), `rotation` = vecteurs
singuliers droits V, `scores` = X_centré V.

## Usage

``` r
pca_fit(X, center = TRUE, scale = FALSE)
```

## Arguments

- X:

  matrice n x p.

- center:

  centrer les colonnes (défaut TRUE).

- scale:

  réduire les colonnes (défaut FALSE).

## Value

liste : `sdev`, `rotation` (p x p), `scores` (n x p), `var_explained`,
`center`, `scale`.
