# Echantillonneur de Metropolis-Hastings (marche aleatoire)

Propose \\x'=x+\mathcal N(0,\text{psd}^2)\\ et l'accepte avec
probabilite \\\min(1,\pi(x')/\pi(x))\\. La chaine converge vers la loi
cible \\\pi\\ (connue a une constante pres). Fonctionne pour toute cible
via sa log-densite.

## Usage

``` r
metropolis_hastings(log_target, init, proposal_sd, n_iter = 10000L)
```

## Arguments

- log_target:

  fonction `x -> log pi(x)` (a une constante additive pres).

- init:

  point de depart (scalaire ou vecteur).

- proposal_sd:

  ecart-type de la proposition (par composante).

- n_iter:

  nombre d'iterations.

## Value

liste : `chain` (n_iter x d), `accept_rate`.
