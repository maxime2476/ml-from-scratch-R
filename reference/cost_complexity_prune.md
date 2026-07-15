# Élagage coût-complexité (éq. 8.4-8.5)

Effondre récursivement les sous-arbres dont le maintien n'améliore pas
\\R\_\alpha(T)=R(T)+\alpha\|T\|\\ (weakest-link, Prop. 8.2).

## Usage

``` r
cost_complexity_prune(object, alpha)
```

## Arguments

- object:

  objet `cart`.

- alpha:

  coût par feuille \\\alpha \ge 0\\.

## Value

objet `cart` élagué.
