# T-learner : CATE par deux modèles séparés

\\\hat\tau(x) = \hat\mu_1(x) - \hat\mu_0(x)\\, chaque \\\mu\\ ajusté sur
le sous-groupe traité / contrôle (forêts du Module 9 par défaut).

## Usage

``` r
t_learner(X, y, d, newX = X, B = 200L)
```

## Arguments

- X:

  data.frame de covariables (apprentissage).

- y:

  réponse.

- d:

  traitement (0/1).

- newX:

  covariables où prédire le CATE (défaut : X).

- B:

  nombre d'arbres.

## Value

vecteur des CATE estimés en newX.
