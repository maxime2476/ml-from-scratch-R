# Sélection des hyperparamètres par maximum de vraisemblance marginale

Maximise la log-vraisemblance marginale (éq. 27.3) sur
\\(\ell,\sigma_f,\sigma_n)\\ — le rasoir d'Occam bayésien automatique.
Optimisation sur l'échelle log (positivité) via `optim` (L-BFGS-B).

## Usage

``` r
gp_optimize(X, y, init = c(1, 1, 0.1))
```

## Arguments

- X, y:

  données.

- init:

  valeurs initiales (échelle, signal, bruit).

## Value

objet `gp` optimisé (avec les hyperparamètres retenus).
