# Tests — Module 13 (théorie de l'apprentissage). Vérifications numériques des
# bornes et objets combinatoires (pas de référence externe : module théorique).

skip_if_not_installed("quadprog")

test_that("hoeffding_bound : formule, bornée par 1, décroissante", {
  expect_equal(hoeffding_bound(100, 0.1), 2 * exp(-2 * 100 * 0.01))
  expect_equal(hoeffding_bound(10, 0), 1)                 # 2*exp(0)=2, plafonné à 1
  expect_lt(hoeffding_bound(200, 0.1), hoeffding_bound(100, 0.1))  # décroît en n
  expect_lt(hoeffding_bound(100, 0.2), hoeffding_bound(100, 0.1))  # décroît en eps
})

test_that("Rademacher empirique linéaire <= borne B*rho/sqrt(n) (éq. 13.5)", {
  set.seed(1)
  for (n in c(40, 150, 600)) {
    X <- matrix(rnorm(n * 4), n, 4)
    est <- empirical_rademacher_linear(X, B = 1, n_draws = 3000, seed = 1)
    expect_lte(est, rademacher_linear_bound(X, B = 1))
  }
})

test_that("Rademacher empirique décroît en ~1/sqrt(n)", {
  set.seed(2)
  X1 <- matrix(rnorm(100 * 4), 100, 4)
  X2 <- matrix(rnorm(400 * 4), 400, 4)
  r1 <- empirical_rademacher_linear(X1, n_draws = 4000, seed = 3)
  r2 <- empirical_rademacher_linear(X2, n_draws = 4000, seed = 3)
  # n x4 -> Rademacher ~ /2 ; on tolère 20 %
  expect_lt(abs(r2 / r1 - 0.5), 0.2)
})

test_that("is_separable : étiquetages séparables et non séparables (XOR)", {
  P4 <- matrix(c(0,0, 1,0, 0,1, 1,1), 4, 2, byrow = TRUE)
  expect_true(is_separable(P4, c(1, 1, -1, -1)))          # demi-plan
  expect_false(is_separable(P4, c(1, -1, -1, 1)))         # XOR : non séparable
})

test_that("shatters_hyperplane : 3 points en position générale oui, 4 non", {
  P3 <- matrix(c(0,0, 1,0, 0,1), 3, 2, byrow = TRUE)
  P4 <- matrix(c(0,0, 1,0, 0,1, 1,1), 4, 2, byrow = TRUE)
  expect_true(shatters_hyperplane(P3))                    # VC >= 3
  expect_false(shatters_hyperplane(P4))                   # VC < 4  (VC = 3)
})

test_that("shatters_hyperplane : 3 points alignés NON pulvérisés", {
  Pc <- matrix(c(0,0, 1,1, 2,2), 3, 2, byrow = TRUE)      # colinéaires
  expect_false(shatters_hyperplane(Pc))                   # (+,-,+) impossible
})

test_that("Hoeffding tient empiriquement (fréquence des déviations <= borne)", {
  set.seed(7)
  n <- 100; p <- 0.5; eps <- 0.15; R <- 20000
  Rhat <- colMeans(matrix(rbinom(n * R, 1, p), n, R))
  freq <- mean(abs(Rhat - p) >= eps)
  expect_lte(freq, hoeffding_bound(n, eps))               # borne respectée
})
