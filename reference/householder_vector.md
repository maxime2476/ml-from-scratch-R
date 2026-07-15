# Vecteur de Householder (choix de signe stable)

Renvoie le vecteur \\v\\ définissant la réflexion \\H(v)=I-2vv^T/v^Tv\\
telle que \\H(v)x = \alpha e_1\\. Implémente l'éq. (0.5) : le signe
\\v_1 = x_1 + \mathrm{sign}(x_1)\\x\\\\ évite l'annulation
catastrophique.

## Usage

``` r
householder_vector(x)
```

## Arguments

- x:

  vecteur non nul.

## Value

le vecteur v (même longueur que x) ; le vecteur nul si x est nul.
