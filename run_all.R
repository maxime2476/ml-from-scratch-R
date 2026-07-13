#!/usr/bin/env Rscript
# =============================================================================
# Orchestrateur de reproduction — ML from Scratch in R
# Reproduit l'ensemble du projet en une commande :
#   Rscript run_all.R [tests|sims|all]
#   - tests : charge R/ et lance toute la suite testthat
#   - sims  : exécute toutes les études Monte Carlo (simulations/)
#   - all   : les deux (défaut)
# Les dérivations et le rapport (.qmd) se rendent séparément via `quarto render`
# ou le livre Quarto (`quarto render` à la racine, cf. _quarto.yml).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args) == 0) "all" else args[1]
root <- normalizePath(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1])), winslash = "/", mustWork = FALSE)
if (is.na(root) || root == "") root <- getwd()
setwd(root)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

run_tests <- function() {
  message("== Chargement des implémentations R/ ==")
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
  suppressMessages(library(testthat))
  if (requireNamespace("mclust", quietly = TRUE)) suppressMessages(library(mclust))
  message("== Suite de tests testthat ==")
  res <- as.data.frame(test_dir("tests/testthat", reporter = "progress",
                                stop_on_failure = FALSE))
  message(sprintf(">> %d passés, %d échecs, %d skip", sum(res$passed),
                  sum(res$failed) + sum(res$error %||% 0), sum(res$skipped)))
  invisible(res)
}

run_sims <- function() {
  sims <- list.files("simulations", pattern = "\\.R$", full.names = TRUE)
  for (s in sims) {
    message("== Simulation : ", basename(s), " ==")
    tryCatch(sys.source(s, envir = new.env()),
             error = function(e) message("  [ERREUR] ", conditionMessage(e)))
  }
}

t0 <- Sys.time()
if (what %in% c("tests", "all")) run_tests()
if (what %in% c("sims", "all"))  run_sims()
message(sprintf("== Terminé en %.1f min ==", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
