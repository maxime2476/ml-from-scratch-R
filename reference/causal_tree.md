# Arbre causal minimal (honnête)

Splits maximisant l'hétérogénéité du traitement \\n_L n_R(\hat\tau_L-
\hat\tau_R)^2\\. Honnêteté : la structure est apprise sur une moitié,
les effets de feuille estimés sur l'autre. Version SIMPLIFIÉE (cf.
dérivation §16.5).

## Usage

``` r
causal_tree(X, y, d, max_depth = 3L, min_leaf = 10L, seed = NULL)
```

## Arguments

- X:

  data.frame de covariables.

- y:

  réponse.

- d:

  traitement (0/1).

- max_depth:

  profondeur maximale.

- min_leaf:

  effectif minimal (par groupe de traitement) par feuille.

- seed:

  graine (partage honnête).

## Value

objet `causal_tree`.
