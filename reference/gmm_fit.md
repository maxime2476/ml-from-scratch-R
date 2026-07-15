# GMM générique (moments non linéaires) par optimisation

Minimise \\\bar g(\theta)^\top W\\\bar g(\theta)\\ (éq. 18.1). Deux
étapes pour la GMM efficace (\\W=\hat S^{-1}\\).

## Usage

``` r
gmm_fit(g_fn, theta0, data, W = NULL, twostep = TRUE)
```

## Arguments

- g_fn:

  fonction `(theta, data) -> matrice n x m` des contributions de moment.

- theta0:

  valeur initiale.

- data:

  données passées à `g_fn`.

- W:

  matrice de pondération initiale (défaut identité).

- twostep:

  GMM efficace à deux étapes (défaut TRUE).

## Value

liste : `coefficients`, `vcov`, `se`, `J`, `J_df`, `J_pvalue`.
