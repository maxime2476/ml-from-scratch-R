# Modele de selection de Heckman (estimation en deux etapes)

Corrige le **biais de selection** : \\y\\ n'est observe que si \\d=1\\
(equation de selection). Etape 1 : probit de \\d\\ sur \\Z_s\\, ratio de
Mills inverse \\\hat\lambda=\phi(Z_s\hat\gamma)/\Phi(Z_s\hat\gamma)\\.
Etape 2 : OLS de \\y\\ sur \\\[X_o,\hat\lambda\]\\ sur les observations
selectionnees. Le coefficient de \\\hat\lambda\\ vaut \\\rho\sigma\\
(nul \\\iff\\ pas de biais de selection).

## Usage

``` r
heckman(selection, outcome, data)
```

## Arguments

- selection:

  formule de selection (LHS binaire).

- outcome:

  formule de resultat (LHS avec NA hors selection).

- data:

  data.frame.

## Value

liste : `gamma` (selection), `beta` (resultat, dont `imr`), `se`.
