# Caractéristiques de Fourier aléatoires (approximation du noyau gaussien)

\\\varphi(x)=\sqrt{2/D}\\\cos(Wx+b)\\, avec \\W\sim\mathcal N(0,\gamma
I)\\ et \\b\sim\mathcal U(0,2\pi)\\. La dimension \\D\\ des
caractéristiques est le **paramètre de complexité** dont on fera varier
la valeur autour de \\n\\.

## Usage

``` r
random_features(X, D, gamma = 1, seed = 1)
```

## Arguments

- X:

  matrice n x p.

- D:

  nombre de caractéristiques.

- gamma:

  échelle (largeur inverse du noyau).

- seed:

  graine (les poids aléatoires doivent être FIXES entre train/test).

## Value

matrice n x D.
