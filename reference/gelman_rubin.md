# Diagnostic de convergence de Gelman-Rubin (statistique R-hat)

Compare la variance INTER-chaines a la variance INTRA-chaine : \\\hat
R=\sqrt{\hat V/W}\\, \\\hat V=\frac{n-1}n W+\frac Bn\\. \\\hat R\to 1\\
a convergence ; \\\>1.1\\ signale une non-convergence.

## Usage

``` r
gelman_rubin(chains)
```

## Arguments

- chains:

  liste de chaines (vecteurs de meme longueur).

## Value

la statistique R-hat.
