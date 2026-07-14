# Tests — Module 40 (classifieurs generatifs). References : MASS, e1071.

make_gen <- function(n = 300, seed = 1) {
  set.seed(seed); p <- 3; K <- 3; mu <- list(c(0, 0, 0), c(2, 2, 0), c(-1, 2, 2))
  X <- do.call(rbind, lapply(1:K, function(k) matrix(rnorm(n / K * p), n / K, p) + rep(mu[[k]], each = n / K)))
  list(X = X, y = factor(rep(1:K, each = n / K)), Xte = matrix(rnorm(200 * p), 200, p) + 1.5)
}

test_that("LDA = MASS::lda (predictions identiques)", {
  skip_if_not_installed("MASS")
  d <- make_gen()
  ph <- lda_predict(lda_fit(d$X, d$y), d$Xte)
  pe <- as.character(predict(MASS::lda(d$X, d$y), d$Xte)$class)
  expect_equal(mean(ph == pe), 1)
})

test_that("QDA = MASS::qda (predictions identiques)", {
  skip_if_not_installed("MASS")
  d <- make_gen()
  ph <- qda_predict(qda_fit(d$X, d$y), d$Xte)
  pe <- as.character(predict(MASS::qda(d$X, d$y), d$Xte)$class)
  expect_equal(mean(ph == pe), 1)
})

test_that("Naive Bayes = e1071::naiveBayes (predictions identiques)", {
  skip_if_not_installed("e1071")
  d <- make_gen()
  ph <- naive_bayes_predict(naive_bayes_fit(d$X, d$y), d$Xte)
  pe <- as.character(predict(e1071::naiveBayes(d$X, d$y), d$Xte))
  expect_equal(mean(ph == pe), 1)
})

test_that("QDA bat LDA quand les covariances different fortement entre classes", {
  set.seed(2); n <- 400
  X1 <- matrix(rnorm(n * 2), n, 2) %*% diag(c(0.3, 3))            # classe 1 : etiree
  X2 <- matrix(rnorm(n * 2), n, 2) %*% diag(c(3, 0.3)) + 1        # classe 2 : etiree l'autre sens
  X <- rbind(X1, X2); y <- factor(rep(1:2, each = n))
  tr <- sample(2 * n, n); te <- setdiff(seq_len(2 * n), tr)
  acc_lda <- mean(lda_predict(lda_fit(X[tr, ], y[tr]), X[te, ]) == y[te])
  acc_qda <- mean(qda_predict(qda_fit(X[tr, ], y[tr]), X[te, ]) == y[te])
  expect_gt(acc_qda, acc_lda)                                     # QDA capte les covariances heterogenes
})

test_that("frontiere LDA lineaire ; QDA quadratique (coherence des discriminants)", {
  d <- make_gen(150)
  # sur l'apprentissage, les trois classent raisonnablement
  expect_gt(mean(lda_predict(lda_fit(d$X, d$y), d$X) == d$y), 0.7)
  expect_gt(mean(qda_predict(qda_fit(d$X, d$y), d$X) == d$y), 0.7)
})
