# Meilleur split d'un nœud (éq. 8.3)

Balayage incrémental sur chaque variable ; renvoie le couple (variable,
seuil) minimisant l'impureté pondérée des enfants (= maximisant \\\Delta
I\\).

## Usage

``` r
best_split(
  X,
  y,
  method,
  kind = "gini",
  min_leaf = 1L,
  classes = NULL,
  mtry = NULL
)
```

## Arguments

- X:

  matrice de prédicteurs (numériques) du nœud.

- y:

  réponse du nœud.

- method:

  "class" ou "anova".

- kind:

  pour la classification : "gini" ou "entropy".

- min_leaf:

  effectif minimal par feuille.

- classes:

  niveaux de classe (classification).

- mtry:

  si non NULL, nombre de variables candidates tirées au hasard pour ce
  split (forêt aléatoire, Module 9) ; NULL = toutes les variables
  (CART).

## Value

liste `gain`, `var` (indice), `val` (seuil) ; `gain = -Inf` si aucun
split.
