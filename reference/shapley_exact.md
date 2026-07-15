# Valeurs de Shapley exactes par énumération (éq. 15.4, p \<= 10)

Fonction de valeur **interventionnelle** \\v(S)=\mathbb
E\_{X_C}\[f(x_S,X_C)\]\\ estimée sur `X_ref`. Complexité \\O(2^p)\\.

## Usage

``` r
shapley_exact(predict_fn, x, X_ref)
```

## Arguments

- predict_fn:

  fonction `data.frame -> numeric`.

- x:

  data.frame d'une ligne (le point à expliquer).

- X_ref:

  data.frame de référence (distribution des variables absentes).

## Value

vecteur nommé des valeurs SHAP (somme `= f(x) - E[f(X_ref)]`).
