# Fonction d'influence d'un MLE (GLM canonique)

\\\mathrm{IC}\_i = n\\\mathcal I^{-1} x_i (y_i-\hat\mu_i)\\, avec
\\\mathcal I = X^\top W X\\ l'information de Fisher (Module 3). Sa
variance \\\frac1{n^2}\sum \mathrm{IC}\_i\mathrm{IC}\_i^\top\\ est le
sandwich robuste ; sous spécification correcte elle vaut \\\mathcal
I^{-1}\\ (variance modèle).

## Usage

``` r
influence_mle(fit)
```

## Arguments

- fit:

  objet `glm_irls` (Module 3, lien canonique).

## Value

liste : `ic` (n x p), `vcov`.
