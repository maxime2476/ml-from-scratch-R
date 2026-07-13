# Laboratoire interactif (Shiny)

Application Shiny pour explorer trois idées centrales du projet, branchées
directement sur les implémentations **from scratch** (`R/`).

## Lancement

```r
# depuis la racine du package
shiny::runApp("shiny-app")
```

## Onglets

1. **Biais-variance (ridge)** — décomposition exacte de l'EQM du ridge
   (`ridge_bias_var`, Module 4, éq. 4.4) le long de la pénalité λ. On visualise
   l'arbitrage biais²/variance et le λ optimal, en faisant varier `n`, `p` et le
   rapport signal/bruit.
2. **Chemins de régularisation** — chemins des coefficients de lasso/ridge
   (`lasso_fit`, `ridge_fit`, Module 4). Le lasso met des coefficients
   exactement à zéro (sélection), le ridge rétrécit sans annuler.
3. **Orthogonalisation (DML)** — modèle partiellement linéaire ; comparaison de
   l'estimateur naïf (biaisé par la confusion) à l'estimateur orthogonal de
   Neyman par partialling-out lasso (Modules 1/4/16). Curseur de force de
   confusion pour voir le biais naïf croître pendant que l'orthogonal reste
   centré sur θ = 1.

Dépendance : `shiny` (dans `Suggests`).
