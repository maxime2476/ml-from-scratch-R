# Prise en main : trois idées, trois lignes

Ce projet réimplémente en R base, depuis leurs dérivations, les modèles
d’apprentissage et d’économétrie. On en donne ici trois aperçus. (En
développement, on charge les sources ; une fois le package installé,
[`library(mlfromscratch)`](https://maxime2476.github.io/ml-from-scratch-R)
suffit.)

``` r

library(mlfromscratch)
```

## 1. Les moindres carrés, à la main = `lm`

``` r

set.seed(1); d <- data.frame(x = rnorm(50)); d$y <- 1 + 2 * d$x + rnorm(50)
maison <- ols_fit(y ~ x, d)$coefficients
reference <- coef(lm(y ~ x, d))
rbind(maison, reference)                       # identiques à la précision machine
#>           (Intercept)        x
#> maison       1.121902 1.954451
#> reference    1.121902 1.954451
```

## 2. Une fonction d’influence *est* un sandwich (Module 24)

La variance « plug-in » de la fonction d’influence de l’OLS **égale**
l’estimateur sandwich HC0 — c’est la même idée sous deux noms.

``` r

X <- model.matrix(y ~ x, d)
inf <- influence_ols(X, d$y)                   # IC de l'OLS
c(influence = sqrt(inf$vcov[2, 2]),
  sandwich  = sqrt(vcov_hc(ols_fit(y ~ x, d), "HC0")[2, 2]))
#> influence  sandwich 
#> 0.1673232 0.1673232
```

## 3. La rétropropagation, par différentiation automatique (Module 28)

Le gradient d’une fonction quelconque, sans dérivée codée à la main :

``` r

ad_grad(function(v) sum(exp(v) * sin(v)), c(0.5, -1.2, 2.0))
#> [1]  2.2373281 -0.1715847  3.6439174
```

## Pour aller plus loin

Chaque module offre une **dérivation** (`derivations/`), une
**implémentation** (`R/`), des **tests** de conformité, et une **étude
Monte Carlo** (`simulations/`). Le fil conducteur — projection, fonction
d’influence, biais-variance, optimisation — est développé dans le
mémoire et la synthèse (`rapport/`).
