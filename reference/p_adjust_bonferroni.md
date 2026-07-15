# Correction de Bonferroni (controle du FWER)

\\\tilde p_i=\min(1,\\m\\p_i)\\. Controle le **taux d'erreur familial**
(probabilite d'au moins un faux positif) au niveau \\\alpha\\, mais est
**conservateur** (perte de puissance quand \\m\\ est grand).

## Usage

``` r
p_adjust_bonferroni(p)
```

## Arguments

- p:

  vecteur de p-valeurs.

## Value

p-valeurs ajustees.
