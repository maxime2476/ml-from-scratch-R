# Tests — Module 39 (SVM). Reference : e1071 (libsvm).

make_svm <- function(n = 120, seed = 1) {
  set.seed(seed); X <- matrix(rnorm(n * 2), n, 2)
  y <- ifelse(X[, 1] + X[, 2] + 0.3 * rnorm(n) > 0, 1, -1)
  list(X = X, y = y)
}

test_that("SVM RBF : predictions et nb de vecteurs de support = e1071", {
  skip_if_not_installed("e1071")
  d <- make_svm(); Xte <- matrix(rnorm(400 * 2), 400, 2)
  m <- svm_fit(d$X, d$y, C = 1, kernel = "rbf", gamma = 0.5)
  e <- e1071::svm(d$X, as.factor(d$y), kernel = "radial", cost = 1, gamma = 0.5, scale = FALSE)
  ph <- svm_predict(m, Xte); pe <- as.numeric(as.character(predict(e, Xte)))
  expect_gt(mean(ph == pe), 0.98)                          # accord quasi parfait
  expect_equal(m$n_sv, e$tot.nSV)                          # meme nombre de SV
})

test_that("SVM lineaire : forte concordance avec e1071", {
  skip_if_not_installed("e1071")
  d <- make_svm(); Xte <- matrix(rnorm(400 * 2), 400, 2)
  m <- svm_fit(d$X, d$y, C = 1, kernel = "linear")
  e <- e1071::svm(d$X, as.factor(d$y), kernel = "linear", cost = 1, scale = FALSE)
  ph <- svm_predict(m, Xte); pe <- as.numeric(as.character(predict(e, Xte)))
  expect_gt(mean(ph == pe), 0.93)
})

test_that("conditions KKT : les alpha respectent 0 <= alpha <= C et sum(alpha y)=0", {
  d <- make_svm(); m <- svm_fit(d$X, d$y, C = 2)
  expect_true(all(m$alpha >= -1e-8 & m$alpha <= 2 + 0.05))
  expect_lt(abs(sum(m$alpha * d$y)), 1e-4)                 # contrainte d'egalite
})

test_that("donnees separables : la marge dure classe parfaitement l'apprentissage", {
  set.seed(2); n <- 80; X <- matrix(rnorm(n * 2), n, 2)
  y <- ifelse(X[, 1] > 0.5, 1, -1); X[y == 1, 1] <- X[y == 1, 1] + 1.5   # marge nette
  m <- svm_fit(X, y, C = 100, kernel = "linear")
  expect_equal(mean(svm_predict(m, X) == y), 1)           # 0 erreur sur l'apprentissage
})

test_that("l'astuce du noyau : RBF separe ce que le lineaire ne peut pas (XOR)", {
  set.seed(3); n <- 200; X <- matrix(runif(n * 2, -1, 1), n, 2)
  y <- ifelse(X[, 1] * X[, 2] > 0, 1, -1)                  # XOR : non lineairement separable
  acc_lin <- mean(svm_predict(svm_fit(X, y, C = 1, "linear"), X) == y)
  acc_rbf <- mean(svm_predict(svm_fit(X, y, C = 1, "rbf", gamma = 2), X) == y)
  expect_lt(acc_lin, 0.8)                                  # le lineaire echoue (~0.5)
  expect_gt(acc_rbf, 0.9)                                  # le RBF reussit
})
