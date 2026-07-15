# Methode delta : variance asymptotique d'une fonction d'un estimateur

Si \\\sqrt n(\hat\theta-\theta)\to\mathcal N(0,\Sigma)\\, alors pour une
fonction reguliere \\g\\, \\\sqrt n(g(\hat\theta)-g(\theta))\to\mathcal
N (0,\nabla g^\top\Sigma\nabla g)\\ (developpement de Taylor au premier
ordre). Le gradient est calcule par differences finies centrees.

## Usage

``` r
delta_method(theta, vcov, g, level = 0.95)
```

## Arguments

- theta:

  estimateur \\\hat\theta\\ (vecteur).

- vcov:

  matrice de covariance estimee \\\hat\Sigma\\.

- g:

  fonction `theta -> scalaire`.

- level:

  niveau de confiance.

## Value

liste : `estimate`, `se`, `ci`.
