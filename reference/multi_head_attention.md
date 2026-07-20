# Attention multi-tetes

Projette \\X\\ en \\Q,K,V\\, decoupe en `n_heads` sous-espaces, applique
l'attention en parallele dans chacun, concatene, puis reprojette. Chaque
tete capte un type de relation different.

## Usage

``` r
multi_head_attention(X, Wq, Wk, Wv, Wo, n_heads = 1L, mask = FALSE)
```

## Arguments

- X:

  entree (T x d_model).

- Wq, Wk, Wv, Wo:

  matrices de projection (d_model x d_model).

- n_heads:

  nombre de tetes (divise d_model)

- mask:

  masque causal.

## Value

liste : `out` (T x d_model), `weights` (liste par tete).
