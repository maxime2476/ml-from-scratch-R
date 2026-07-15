# Estimateur 2SLS (forme fermée, éq. 5.4)

Projette X sur l'espace des instruments (\\\hat X = P_Z X\\, étape 1)
puis régresse y sur \\\hat X\\ (étape 2), d'où \\\hat\beta = (X^T P_Z
X)^{-1} X^T P_Z y\\. La variance (éq. 5.5) utilise les résidus sur le
VRAI X (équation structurelle), pas \\\hat X\\.

## Usage

``` r
tsls_fit(y, X, Z)
```

## Arguments

- y:

  réponse (vecteur longueur n).

- X:

  matrice n x k des régresseurs (constante incluse ; colonnes exogènes
  et endogènes).

- Z:

  matrice n x m des instruments (m \>= k ; inclut la constante et les
  régresseurs exogènes, plus les instruments exclus).

## Value

liste : `coefficients`, `vcov`, `se`, `sigma2`, `df.residual`, `fitted`,
`residuals`, `Xhat`.
