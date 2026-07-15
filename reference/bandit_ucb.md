# Bandit UCB1 (borne de confiance superieure)

Choisit a chaque tour le bras maximisant \\\hat\mu_a+\sqrt{2\log
t/n_a}\\ : l'optimisme face a l'incertitude. Regret **logarithmique**
\\O(\log T)\\.

## Usage

``` r
bandit_ucb(means, horizon, seed = NULL)
```

## Arguments

- means:

  moyennes VRAIES des bras (Bernoulli ou bornees) ; @param horizon T ;

- seed:

  graine.

## Value

liste : `regret` (cumule), `arms` (choisis), `counts`.
