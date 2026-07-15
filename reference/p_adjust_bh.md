# Correction de Benjamini-Hochberg (controle du FDR)

Controle le **taux de fausses decouvertes** (proportion attendue de faux
positifs PARMI les rejets) au niveau \\\alpha\\ : trie les p-valeurs,
applique \\\tilde p\_{(i)}=\min\_{j\ge i} m\\p\_{(j)}/j\\. Bien plus
**puissant** que Bonferroni quand beaucoup d'hypotheses sont vraiment
fausses.

## Usage

``` r
p_adjust_bh(p)
```

## Arguments

- p:

  vecteur de p-valeurs.

## Value

p-valeurs ajustees.
