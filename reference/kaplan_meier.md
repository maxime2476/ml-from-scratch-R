# Estimateur de survie de Kaplan-Meier

\\\hat S(t)=\prod\_{t_i\le t}\bigl(1-d_i/n_i\bigr)\\, ou \\d_i\\ deces
et \\n_i\\ sujets a risque en \\t_i\\. Gere la **censure** a droite (les
censures reduisent le risque sans compter comme evenement).

## Usage

``` r
kaplan_meier(time, event)
```

## Arguments

- time:

  durees observees

- event:

  1 = evenement, 0 = censure.

## Value

liste : `time` (temps d'evenement), `surv`, `n_risk`, `n_event`.
