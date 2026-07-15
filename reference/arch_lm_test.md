# Test ARCH-LM (heteroscedasticite conditionnelle)

Regresse \\\hat\varepsilon_t^2\\ sur ses \\q\\ retards ; \\nR^2\sim
\chi^2_q\\ sous absence d'effet ARCH. Rejet = la volatilite se regroupe.

## Usage

``` r
arch_lm_test(x, q = 5L)
```

## Arguments

- x:

  serie (rendements) ; @param q nombre de retards.

## Value

liste : `statistic`, `df`, `p_value`.
