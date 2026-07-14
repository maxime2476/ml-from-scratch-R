# =============================================================================
# Runner de tests robuste — un PROCESSUS FRAIS par fichier de test.
# -----------------------------------------------------------------------------
# La suite valide chaque module contre son package de référence (glmnet, AER,
# gmm, DoubleML, grf, iml, hdm, sensemakr...). Chargés ENSEMBLE dans un même
# processus R, ces packages lourds se marchent dessus (S4 setGeneric sur
# `coef`/`mean`, symboles exportés, état C accumulé), ce qui corrompt
# l'environnement et fait échouer des tests SANS RAPPORT avec le code testé —
# alors que chaque fichier passe seul. Ce runner exécute donc chaque fichier de
# test dans un sous-processus R FRAIS (via callr), garantissant l'isolation.
#
# Usage :  Rscript run_tests.R
# =============================================================================

stopifnot(requireNamespace("callr", quietly = TRUE),
          requireNamespace("testthat", quietly = TRUE))
root <- normalizePath(".")
files <- sort(list.files(file.path(root, "tests", "testthat"),
                         pattern = "^test-.*[.]R$", full.names = TRUE))

run_one <- function(f, root) {
  callr::r(function(file, root) {
    setwd(file.path(root, "tests"))
    library(testthat)
    res <- as.data.frame(test_file(file, reporter = "silent", stop_on_failure = FALSE))
    c(pass = sum(res$passed), fail = sum(res$failed), error = sum(res$error))
  }, args = list(file = f, root = root), show = FALSE)
}

cat(sprintf("Exécution de %d fichiers de test en processus isolés...\n\n", length(files)))
tot <- c(pass = 0, fail = 0, error = 0)
for (f in files) {
  r <- tryCatch(run_one(f, root),
                error = function(e) { cat("  [CRASH]", basename(f), ":", conditionMessage(e), "\n"); c(pass = 0, fail = 0, error = 1) })
  tot <- tot + r
  flag <- if (r["fail"] + r["error"] > 0) "  <-- ÉCHEC" else ""
  cat(sprintf("  %-22s pass %3d  fail %d  error %d%s\n", basename(f), r["pass"], r["fail"], r["error"], flag))
}
cat(sprintf("\n===== TOTAL : pass %d | fail %d | error %d =====\n", tot["pass"], tot["fail"], tot["error"]))
if (tot["fail"] + tot["error"] > 0) quit(status = 1)
