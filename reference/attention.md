# Attention par produit scalaire mis a l'echelle

\\\mathrm{Attn}(Q,K,V)=\mathrm{softmax}\\\bigl(QK^\top/\sqrt{d_k}\bigr)V\\.
Chaque requete \\q_i\\ interroge toutes les cles ; les poids
(similarites normalisees) melangent les valeurs. Une option de masque
**causal** interdit de regarder le futur (auto-regression).

## Usage

``` r
attention(Q, K, V, mask = FALSE)
```

## Arguments

- Q, K, V:

  matrices requetes (T_q x d_k), cles (T_k x d_k), valeurs (T_k x d_v).

- mask:

  logique : masque causal (position i ne voit que j \<= i).

## Value

liste : `out` (T_q x d_v), `weights` (T_q x T_k).
