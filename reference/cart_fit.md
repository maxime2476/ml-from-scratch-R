# Ajustement d'un arbre CART

Croissance récursive gloutonne (§8.3-8.4) : à chaque nœud, `best_split`
choisit le split optimal ; arrêt sur profondeur, effectif, pureté ou
gain nul.

## Usage

``` r
cart_fit(
  formula,
  data,
  method = c("class", "anova"),
  kind = "gini",
  max_depth = 30L,
  min_split = 20L,
  min_leaf = 7L,
  min_gain = 1e-09,
  mtry = NULL
)
```

## Arguments

- formula:

  formule façon `rpart` (prédicteurs numériques).

- data:

  data.frame.

- method:

  "class" (Gini/entropie) ou "anova" (variance).

- kind:

  impureté de classification : "gini" (défaut) ou "entropy".

- max_depth:

  profondeur maximale.

- min_split:

  effectif minimal pour tenter un split.

- min_leaf:

  effectif minimal par feuille.

- min_gain:

  gain d'impureté minimal pour accepter un split.

- mtry:

  variables candidates par split (forêt aléatoire) ; NULL = toutes.

## Value

objet `cart` (arbre + métadonnées).
