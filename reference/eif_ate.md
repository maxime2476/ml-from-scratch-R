# Fonction d'influence EFFICACE de l'ATE (doublement robuste / AIPW)

\\\psi_i = \mu_1(x_i)-\mu_0(x_i) + \frac{d_i(y_i-\mu_1)}{e_i} -
\frac{(1-d_i)(y_i-\mu_0)}{1-e_i} - \tau\\. Son second moment est la
**borne d'efficacité semiparamétrique** de l'ATE ; l'AIPW/DML (Module
16) l'atteint.

## Usage

``` r
eif_ate(y, d, mu1, mu0, e)
```

## Arguments

- y:

  résultat.

- d:

  traitement (0/1).

- mu1, mu0:

  espérances conditionnelles estimées \\E\[Y\|X,D=1/0\]\\.

- e:

  score de propension estimé \\P(D=1\|X)\\.

## Value

liste : `ate`, `eif`, `se` (= sqrt(borne d'efficacité / n)).
