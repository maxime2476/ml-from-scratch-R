# Agrégation des ATT(g,t) (Callaway-Sant'Anna)

Agrégation des ATT(g,t) (Callaway-Sant'Anna)

## Usage

``` r
aggregate_att(attgt, type = c("simple", "dynamic", "group"))
```

## Arguments

- attgt:

  sortie de `att_gt`.

- type:

  "simple" (ATT global, périodes post, pondéré par taille de groupe),
  "dynamic" (event-study par ancienneté \\e=t-g\\) ou "group" (par
  cohorte).

## Value

pour "simple" un scalaire ; sinon un data.frame (`e`/`g`, `att`).
