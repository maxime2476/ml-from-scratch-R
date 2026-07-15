# Selection de fenetre par validation croisee leave-one-out (Nadaraya-Watson)

Minimise \\\sum_i (y_i-\hat m\_{-i}(x_i))^2\\, ou \\\hat m\_{-i}\\
exclut l'observation \\i\\ (formule du residu : retirer le poids
diagonal).

## Usage

``` r
bw_loocv(x, y, bws)
```

## Arguments

- x, y:

  donnees

- bws:

  fenetres candidates

## Value

liste : `bw` (optimale), `cv` (erreurs CV par fenetre).
