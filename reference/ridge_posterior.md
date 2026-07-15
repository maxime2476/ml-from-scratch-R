# Postérieure conjuguée du ridge (Prop. 14.2)

Sous \\y\mid\beta\sim\mathcal N(X\beta,\sigma^2 I)\\ et prior
\\\beta\sim\mathcal N(0,\tau^2 I)\\, renvoie la moyenne et la covariance
a posteriori (éq. 14.3), avec \\\lambda=\sigma^2/\tau^2\\.

## Usage

``` r
ridge_posterior(X, y, lambda, sigma2)
```

## Arguments

- X:

  matrice de design n x p (sans intercept, ou déjà centré).

- y:

  réponse.

- lambda:

  pénalité ridge \\\lambda=\sigma^2/\tau^2\\.

- sigma2:

  variance du bruit (connue).

## Value

liste : `mean` (= estimateur ridge), `cov` (\\\Sigma\_{\text{post}}\\).
