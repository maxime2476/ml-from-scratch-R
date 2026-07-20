# Q-learning tabulaire (MDP a modele INCONNU)

Apprend \\Q^\\\\ par interaction, sans connaitre \\P,R\\ : politique
\\\varepsilon\\-greedy, mise a jour par difference temporelle
\\Q(s,a)\leftarrow Q(s,a)+\alpha\[r+\gamma\max\_{a'}Q(s',a')-Q(s,a)\]\\.
Converge vers l'optimum du `value_iteration`.

## Usage

``` r
q_learning(
  P,
  R,
  gamma = 0.9,
  episodes = 3000L,
  steps = 30L,
  alpha = 0.1,
  epsilon = 0.1,
  seed = NULL
)
```

## Arguments

- P, R, gamma:

  cf. `value_iteration` (utilises pour SIMULER l'environnement) ;

- episodes, steps:

  duree

- alpha, epsilon:

  taux et exploration

- seed:

  graine.

## Value

liste : `Q`, `policy`.
