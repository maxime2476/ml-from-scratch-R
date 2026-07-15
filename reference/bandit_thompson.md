# Bandit par echantillonnage de Thompson (Bernoulli)

Prior Beta(1,1) par bras ; a chaque tour, tire \\\theta_a\sim\text{Beta}
(\alpha_a,\beta_a)\\ et joue \\\arg\max\theta_a\\ ; met a jour la
posterieure. Probability matching bayesien, regret logarithmique.

## Usage

``` r
bandit_thompson(means, horizon, seed = NULL)
```

## Arguments

- means, horizon, seed:

  cf. `bandit_ucb`.

## Value

liste : `regret`, `arms`, `counts`.
