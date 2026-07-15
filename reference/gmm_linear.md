# GMM linéaire (moments d'instruments) : 2SLS et GMM efficace à deux étapes

Moments \\g_i(\beta)=Z_i(y_i-X_i^\top\beta)\\. Avec `twostep = FALSE` et
\\W=(Z^TZ)^{-1}\\, renvoie le 2SLS (Prop. 18.2). Avec `twostep = TRUE`,
renvoie la GMM efficace (pondération \\\hat S^{-1}\\, robuste à
l'hétéroscédasticité) et le test de suridentification J (éq. 18.3).

## Usage

``` r
gmm_linear(y, X, Z, twostep = TRUE)
```

## Arguments

- y:

  réponse.

- X:

  régresseurs (matrice n x k, constante incluse).

- Z:

  instruments (matrice n x m, m \>= k).

- twostep:

  GMM efficace à deux étapes (défaut TRUE) ; sinon 2SLS.

## Value

liste : `coefficients`, `vcov`, `se`, `J`, `J_df`, `J_pvalue`, `step`.
