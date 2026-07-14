# =============================================================================
# Pipeline `targets` — reproduction à dépendances suivies (compendium citable)
# -----------------------------------------------------------------------------
# `targets` construit un graphe de dépendances : seuls les objets dont une
# dépendance a changé sont recalculés (make pour R). Ici :
#   R/*.R  ->  chargement des fonctions
#   simulations/*.R  ->  chaque étude Monte Carlo (branche dynamique)
#   tests  ->  suite isolée (run_tests.R)
# Usage :  targets::tar_make()   |   targets::tar_visnetwork()  (graphe)
# =============================================================================

library(targets)

# fonction exécutant un script de simulation dans un environnement frais
run_sim <- function(path) {
  sys.source(path, envir = new.env(parent = globalenv()))
  path
}

list(
  # --- sources suivies comme fichiers (recalcul si un .R change) ------------
  tar_target(r_files,   list.files("R", pattern = "[.]R$", full.names = TRUE), format = "file"),
  tar_target(sim_files, list.files("simulations", pattern = "[.]R$", full.names = TRUE), format = "file"),

  # --- chargement des implémentations ---------------------------------------
  tar_target(functions, { for (f in r_files) source(f); length(r_files) }),

  # --- suite de tests, un processus par fichier (isolation) -----------------
  tar_target(tests, {
    functions
    system2("Rscript", "run_tests.R", stdout = TRUE, stderr = TRUE)
  }),

  # --- études Monte Carlo : une branche par fichier -------------------------
  tar_target(simulations, run_sim(sim_files),
             pattern = map(sim_files), format = "file"),

  # --- livre Quarto (dérivations + rapport) ---------------------------------
  tar_target(book, {
    simulations                                   # les figures doivent exister
    if (nzchar(Sys.which("quarto"))) system2("quarto", "render") else "quarto absent"
  })
)
