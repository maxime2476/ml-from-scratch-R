# Prediction d'un SVM

\\\hat y(x)=\mathrm{sign}\bigl(\sum_i\alpha_i y_i K(x_i,x)+b\bigr)\\.

## Usage

``` r
svm_predict(object, Xnew, decision = FALSE)
```

## Arguments

- object:

  objet `svm_fit` ; @param Xnew nouvelles observations ;

- decision:

  renvoyer la valeur de decision (defaut FALSE : le signe).

## Value

vecteur de classes (ou de valeurs de decision).
