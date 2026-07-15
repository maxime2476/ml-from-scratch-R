# Prédiction d'un processus gaussien

Moyenne a posteriori \\\bar f\_\*=k\_\*^\top\alpha\\ et variance de la
**fonction latente**
\\\operatorname{Var}\[f\_\*\]=k\_{\*\*}-k\_\*^\top(K+\sigma_n^2I)^{-1}k\_\*\\.
La variance d'une **observation** bruitée ajoute \\\sigma_n^2\\.

## Usage

``` r
gp_predict(object, Xnew)
```

## Arguments

- object:

  objet `gp`.

- Xnew:

  points de prédiction.

## Value

liste : `mean`, `sd` (fonction latente), `sd_obs` (observation bruitée).
