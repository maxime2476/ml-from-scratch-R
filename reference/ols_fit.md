# Ajustement MCO par QR, avec inférence

Résout les équations normales (éq. 1.2) via la QR du Module 0, puis
calcule la matrice de variance classique \\s^2(X^TX)^{-1}\\ (éq. 1.5),
l'estimateur sans biais \\s^2\\ (éq. 1.6) et les leviers (diagonale de
la matrice chapeau). Aucune inversion de \\X^TX\\ n'est formée : on
exploite le facteur \\R\\ de la QR (\\X^TX = R^T R\\, donc
\\(X^TX)^{-1}=R^{-1}R^{-T}\\).

## Usage

``` r
ols_fit(formula, data)
```

## Arguments

- formula:

  formule façon `lm` (ex. `y ~ x1 + x2`).

- data:

  data.frame contenant les variables de la formule.

## Value

objet de classe `ols` : `coefficients`, `vcov`, `sigma2`, `df.residual`,
`fitted`, `residuals`, `hat`, `rss`, `tss`, `rank`, `n`, `p`,
`has_intercept`, `terms`, `xlevels`.
