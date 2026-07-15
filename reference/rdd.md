# Discontinuite de regression (RDD, effet local au seuil)

Ajuste une regression **locale lineaire** de chaque cote d'un seuil
\\c\\ et estime le **saut** \\\tau=\lim\_{x\downarrow
c}m(x)-\lim\_{x\uparrow c}m(x)\\, l'effet causal local (design
quasi-experimental de Thistlethwaite-Campbell).

## Usage

``` r
rdd(y, x, cutoff = 0, bw, kernel = c("tri", "gauss"))
```

## Arguments

- y:

  resultat

- x:

  variable de forcage (running variable)

- cutoff:

  seuil \\c\\

- bw:

  fenetre (des deux cotes)

- kernel:

  "tri" (triangulaire, defaut) ou "gauss".

## Value

liste : `tau` (saut), `mu_left`, `mu_right`, `n_left`, `n_right`.
